
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
      "porothermoelastic.alpha_d",
      "porothermoelastic.young",
      "porothermoelastic.poisson",
      "porothermoelastic.lame_mu",
      "porothermoelastic.lame_lambda"
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
  // Outputs in the ElementCoupler
  auto & dbl_params = trg_ec.dbl_params;
  auto & vector_params = trg_ec.vector_params;
  auto & tensor_params = trg_ec.tensor_params;
  vector<double> & von_mises      = dbl_params["von_mises"];     von_mises.clear();
  vector<double> & epskk          = dbl_params["epsilon"];       epskk.clear();
  vector<RealTensor> & sigeff     = tensor_params["sigeff"];     sigeff.clear(); 
  vector<RealTensor> & sigtot     = tensor_params["sigtot"];     sigtot.clear();
  vector<RealTensor> & deviatoric = tensor_params["deviatoric"]; deviatoric.clear();

  uint nqp = U.size();

  for (uint qp=0; qp<nqp; qp++) 
  {
    double epskk_ = 0;
    for (uint k=0; k<3; k++ ) epskk_ += GRAD_U[qp](k,k);

    RealTensor sigeff_;
    for (uint i=0; i<3; i++) for (uint j=0; j<3; j++) 
    for (uint k=0; k<3; k++) for (uint l=0; l<3; l++)
      sigeff_(i,j) += C_ijkl(qp,i,j,k,l) * GRAD_U[qp](k,l);

    // \sig_tot = sig_eff - \alpha_d * T * \delta_ij
    RealTensor sigtot_ = sigeff_;
    for (uint k=0; k<3; k++ ) sigtot_(k,k) -= alpha_d[qp] * temperature[qp];

    // dev_ij = sig_ij - 1/3 \delta_ij \sigma_kk
    RealTensor deviatoric_ = sigtot_;
    for (uint i=0; i<3; i++ )
    for (uint k=0; k<3; k++ ) 
      deviatoric_(i,i) -= (1./3.) * sigtot_(k,k);

    double J2=0;
    for (uint i=0; i<3; i++ )
    for (uint j=0; j<3; j++ ) 
      J2 += (1./2.) * deviatoric_(i,j) * deviatoric_(i,j);
    double von_mises_ = sqrt( 3 * J2 );

    von_mises.push_back( von_mises_ );
    deviatoric.push_back( deviatoric_ );
    epskk.push_back( epskk_ );
    sigeff.push_back( sigeff_ );
    sigtot.push_back( sigtot_ );
  }
}



/**
 *
 */
void StressMaterial::setup_variables()
{
  SCOPELOG(1);

  // TODO: fetch this info from configuration
  Order order = SECOND;
  FEFamily fef = L2_LAGRANGE;
//  if ( ! order ) fef = MONOMIAL;  // a constant is a monomial

  vector<string> sname = { "sigeff", "sigtot", "deviatoric" };
  vector<string> sdir  = { "XX",  "YY",  "ZZ",  "XY",  "XZ",   "YZ" };

  set<subdomain_id_type> sids = { sid };
  for ( auto sn : sname ) 
  for ( auto sd : sdir )
    system.add_variable(sn+sd, order, fef, &sids);

  system.add_variable("von_mises", order, fef, &sids);
}

/**
 *
 */
void StressMaterial::init_fem()
{
  SCOPELOG(1);
  uint vid = system.variable_number( "sigtotXX" );
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

  get_from_element_coupler(  "porothermoelastic.alpha_d",       alpha_d        ); 
  get_from_element_coupler(  "porothermoelastic.lame_mu",       lame_mu        ); 
  get_from_element_coupler(  "porothermoelastic.lame_lambda",   lame_lambda    ); 

  // From the solvers
  get_from_element_coupler(  "T",             temperature    ); 
  get_from_element_coupler(  "U",             U              ); 
  get_from_element_coupler(  "GRAD_U",        GRAD_U         ); 
}
