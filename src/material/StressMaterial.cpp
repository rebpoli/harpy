
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
  Material( sid_, config ), system( sys_ ), U(3)
{
  SCOPELOG(1);

  // Lists the necessary properties to fetch from the config during init_coupler
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
 *   Do the actual stress calculations
 */
void StressMaterial::feed_coupler( ElemCoupler & trg_ec )
{
  SCOPELOG(1);
  vector<string> vnames = { "sigTOTXX", "sigTOTYY", "sigTOTZZ", "sigTOTXY", "sigTOTXZ", "sigTOTYZ" };

  const vector<vector<RealGradient>> & dphi = fe->get_dphi();
  const vector<vector<Real>> & phi = fe->get_phi();
  const vector<Real> & jxw = fe->get_JxW();
  DofMap & dof_map = system.get_dof_map();
  
  uint nqp = qrule.n_points();

  // Initialize outputs in the ElementCoupler
  auto & dbl_params = trg_ec.dbl_params;
  auto & vector_params = trg_ec.vector_params;
  auto & tensor_params = trg_ec.tensor_params;
  vector<double> & von_mises    = dbl_params["von_mises"]; von_mises.clear();
  vector<double> & epskk        = dbl_params["epsilon"]; epskk.clear();
  vector<RealTensor> & sigeff   = tensor_params["sigeff"]; sigeff.clear();
  vector<RealTensor> & sigtot   = tensor_params["sigtot"]; sigtot.clear();

  for (uint qp=0; qp<nqp; qp++) 
  for (uint k=0; k<3; k++ )
    epskk[qp] += GRAD_U[qp](k,k);

  for (uint qp=0; qp<nqp; qp++) 
  for (uint i=0; i<3; i++) for (uint j=0; j<3; j++) 
  for (uint k=0; k<3; k++) for (uint l=0; l<3; l++)
    sigeff[qp](i,j) += C_ijkl(qp,i,j,k,l) * GRAD_U[qp](k,l);

  const auto kd = Math::kronecker_delta;
  for (uint qp=0; qp<nqp; qp++) 
  for (uint i=0; i<3; i++)
  for (uint j=0; j<3; j++) 
    sigeff[qp](i,j) += sigeff[qp](i,j);
}

/**
 *
 */
void StressMaterial::project( ElemCoupler & ec, string vname )
{
  const std::vector<std::vector<Real>> & phi = fe->get_phi();
  const std::vector<Real> & jxw = fe->get_JxW();
  const vector<double> & vals_qp = ec.dbl_params[vname];
  uint vid = system.variable_number(vname);
  Elem & elem = system.get_mesh().elem_ref( ec.eid );

  //  M x = F
  uint n_dofs = phi.size();
  DenseMatrix<Real> MAT(n_dofs, n_dofs);
  DenseVector<Number> F(n_dofs);

  uint nqp = qrule.n_points();
  for(uint qp=0; qp<nqp; qp++) 
  for(uint B=0; B<n_dofs; B++) 
    F(B) += jxw[qp] * vals_qp[qp] * phi[B][qp];

  for(uint qp=0; qp<nqp; qp++) 
  for(uint B=0; B<n_dofs; B++) 
  for(uint M=0; M<n_dofs; M++) 
    MAT(B,M) += jxw[qp] * phi[B][qp] * phi[M][qp];

  for(uint B=0; B<n_dofs; B++) 
    if ( std::abs(MAT(B,B)) < 1e-10 ) MAT(B,B) = 1;

  DenseVector<Number> X;
//  MAT.cholesky_solve(F, X);
  MAT.lu_solve(F, X);

  const DofMap & dof_map = system.get_dof_map();
  std::vector<dof_id_type> dof_indices;
  dof_map.dof_indices (&elem, dof_indices, vid);

  // Cada processador seta os seus nos apenas
  dof_id_type f = system.solution->first_local_index();
  dof_id_type l = system.solution->last_local_index();
  for(uint B=0; B<n_dofs; B++)
  {
    dof_id_type dof_i = dof_indices[B];
    if ((f <= dof_i) && (dof_i < l )) 
      system.solution->set(dof_i, X(B));
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
 *
 *
 */
void StressMaterial::reinit( const Elem & elem, uint side ) 
{
  fe->reinit( &elem );

  get_from_element_coupler(  "T",             temperature    ); 
  get_from_element_coupler(  "alpha_d",       alpha_d        ); 
  get_from_element_coupler(  "lame_mu",       lame_mu        ); 
  get_from_element_coupler(  "lame_lambda",   lame_lambda    ); 
  get_from_element_coupler(  "U",             U              ); 
  get_from_element_coupler(  "GRAD_U",        GRAD_U         ); 
}
