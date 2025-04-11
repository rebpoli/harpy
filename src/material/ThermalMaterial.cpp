
#include "material/ThermalMaterial.h"
#include "config/MaterialConfig.h"

#include "libmesh/system.h"
#include "libmesh/equation_systems.h"

#include "libmesh/string_to_enum.h"
#include "libmesh/enum_fe_family.h"
#include "libmesh/enum_order.h"

#include "libmesh/explicit_system.h"

#include "util/OutputOperators.h"


/**
 *
 */
ThermalMaterial::ThermalMaterial( suint sid_,
                                  const MaterialConfig & config, 
                                  ExplicitSystem & sys_ ) :
  Material( sid_, config ), system( sys_ )
{
  SCOPELOG(1);
  dlog(1) << config;
  setup_variables();
}

/**
 *
 */
ThermalMaterial::~ThermalMaterial()
{ 
SCOPELOG(1);
}

/**
 *
 */
void ThermalMaterial::project( ElemCoupler & ec, string vname )
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
void ThermalMaterial::setup_variables()
{
  set<subdomain_id_type> sids = { sid };
  uint vid = system.add_variable( "T", SECOND, L2_LAGRANGE, &sids );
  dlog(1) << "Adding variable T - sid=" << sid << "";
}

/**
 *
 */
void ThermalMaterial::init_fem()
{
  SCOPELOG(1);
  uint vid = system.variable_number( "T" );
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
void ThermalMaterial::reinit( const Elem & elem, uint side ) 
{
  fe->reinit( &elem );
}
