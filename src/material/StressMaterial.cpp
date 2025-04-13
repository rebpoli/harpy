
#include "material/StressMaterial.h"

#include "config/MaterialConfig.h"
#include "harpy/Coupler.h"
#include "util/OutputOperators.h"

#include "libmesh/system.h"
#include "libmesh/equation_systems.h"
#include "libmesh/explicit_system.h"

/**
 *
 */
StressMaterial::StressMaterial( suint sid_,
                                  const MaterialConfig & config, 
                                  ExplicitSystem & sys_ ) :
  MaterialExplicit( sid_, config, sys_ ), U(3)
{
  SCOPELOG(1);

  // Lists the necessary properties to fetch from the config during Materia::init_coupler
  required_material_properties.assign({
      "alpha_d", "young", "poisson", "lame_mu", "lame_lambda"
  });

  dlog(1) << config;
  setup_variables();
}

/**
 *
 */
StressMaterial::~StressMaterial()
{ SCOPELOG(1); }

/**
 *     We assume the src and target element couplers are in the same mesh.
 *
 *     In the future we may be interested in interpolating stresses, but not for now.
 */
void StressMaterial::feed_coupler( ElemCoupler & trg_ec )
{
  // Inputs
  auto & U = elem_coupler->vector_params.at("U");
  auto & GRAD_U = elem_coupler->tensor_params.at("GRAD_U");

  // Outputs in the ElementCoupler
  auto & dbl_params = trg_ec.dbl_params;
  auto & vector_params = trg_ec.vector_params;
  auto & tensor_params = trg_ec.tensor_params;
  vector<double> & von_mises    = dbl_params["von_mises"];  von_mises.clear();
  vector<double> & epskk        = dbl_params["epsilon"];  epskk.clear();
  vector<RealTensor> & sigeff   = tensor_params["sigeff"];  sigeff.clear(); 
  vector<RealTensor> & sigtot   = tensor_params["sigtot"]; sigtot.clear();


  uint nqp = U.size();

  for (uint qp=0; qp<nqp; qp++) 
  {
    double epskk_ = 0;
    for (uint k=0; k<3; k++ ) epskk_ += GRAD_U[qp](k,k);

    RealTensor sigeff_;
    for (uint i=0; i<3; i++) for (uint j=0; j<3; j++) 
    for (uint k=0; k<3; k++) for (uint l=0; l<3; l++)
      sigeff_(i,j) += C_ijkl(qp,i,j,k,l) * GRAD_U[qp](k,l);

    epskk.push_back( epskk_ );
    sigeff.push_back( sigeff_ );
  }
}



/**
 *
 */
void StressMaterial::setup_variables()
{
  SCOPELOG(1);

  // TODO: fetch this info from configuration
  Order order = FIRST;
  FEFamily fef = L2_LAGRANGE;
//  if ( ! order ) fef = MONOMIAL;  // a constant is a monomial

  set<subdomain_id_type> sids = { sid };
  system.add_variable("sigTOTXX", order, fef, &sids);
  system.add_variable("sigTOTYY", order, fef, &sids);
  system.add_variable("sigTOTZZ", order, fef, &sids);
  system.add_variable("sigTOTXY", order, fef, &sids);
  system.add_variable("sigTOTXZ", order, fef, &sids);
  system.add_variable("sigTOTYZ", order, fef, &sids);

  system.add_variable("sigeffXX", order, fef, &sids);
  system.add_variable("sigeffYY", order, fef, &sids);
  system.add_variable("sigeffZZ", order, fef, &sids);
  system.add_variable("sigeffXY", order, fef, &sids);
  system.add_variable("sigeffXZ", order, fef, &sids);
  system.add_variable("sigeffYZ", order, fef, &sids);
}

/**
 *
 */
void StressMaterial::init_fem()
{
  SCOPELOG(1);
  uint vid = system.variable_number( "sigTOTXX" );
  DofMap & dof_map = system.get_dof_map();
  FEType fe_type = dof_map.variable_type(vid);

  fe = move( FEBase::build(3, fe_type) );

  qrule = QGauss( 3, fe_type.default_quadrature_order() );

  fe->attach_quadrature_rule (&qrule);

  // Enable calculations calculations
  fe->get_JxW();
  fe->get_phi();
  fe->get_dphi();
  fe->get_xyz();
}

/**
 *    The reinit needs a coupler, that is the source of information
 *    for the stress calculations.
 *
 *    The StressSolver (and StressMaterial) operate between couplers:
 *    
 *    - the VPSolver feeds the StressSolver.coupler
 *    - The Stress material copies the information from the coupler to its object
 *    - The feed_coupler uses this information to feed the VPSolver coupler back with the stresses
 *
 *    NOTE: the couplers are individualized for each solver because the mesh _might_ be different
 */
void StressMaterial::reinit( Coupler & coupler, const Elem & elem, uint side ) 
{
  // Initialize FE and the element coupler (the inputs for the computations)
  fe->reinit( &elem );
  elem_coupler = & ( coupler.elem_coupler( elem.id() ) );

  get_from_element_coupler(  "T",             temperature    ); 
  get_from_element_coupler(  "alpha_d",       alpha_d        ); 
  get_from_element_coupler(  "lame_mu",       lame_mu        ); 
  get_from_element_coupler(  "lame_lambda",   lame_lambda    ); 
  get_from_element_coupler(  "U",             U              ); 
  get_from_element_coupler(  "GRAD_U",        GRAD_U         ); 
}
