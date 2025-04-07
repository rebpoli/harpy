
#include "material/ViscoPlasticMaterial.h"
#include "config/MaterialConfig.h"

#include "libmesh/system.h"
#include "libmesh/equation_systems.h"

#include "libmesh/string_to_enum.h"
#include "libmesh/enum_fe_family.h"
#include "libmesh/enum_order.h"


/**
 *
 */
ViscoPlasticMaterial::ViscoPlasticMaterial( suint sid_, const MaterialConfig & config, System & sys_ ) :
  Material(sid_, config, sys_),
  dof_indices_var(3)
{
  setup_variables();
}

/**
 *
 */
void ViscoPlasticMaterial::setup_variables()
{
  // Add displacements variables to the current subdomain ID
  if (! config.fem_by_var.count( "U" ) ) flog << "Undefined var setup for variable 'U'. Please revise model material.FEM section.";
  auto & femspec = config.fem_by_var.at("U");
  Order order = Utility::string_to_enum<Order>( femspec.order ) ;
  FEFamily fe_family = Utility::string_to_enum<FEFamily>( femspec.family );

  dlog(1) << "Setting up variable 'U' for ViscoPlasticMaterial (sid=" << sid << ") ...";
  dlog(1) << "     Order:" << order;
  dlog(1) << "     FEFamily:" << fe_family;

  set<subdomain_id_type> sids = { sid };
  system.add_variable( "UX", order, fe_family, &sids );
  system.add_variable( "UY", order, fe_family, &sids );
  system.add_variable( "UZ", order, fe_family, &sids );
}

/**
 *  This can only be done after EquationSystems init...
 */
void ViscoPlasticMaterial::init_fem()
{
  SCOPELOG(1);
  uint vid = system.variable_number( "UX" );

  // Setup shape functions
  uint dim = 3;
  DofMap & dof_map = system.get_dof_map();
  FEType fe_type = dof_map.variable_type(vid);
  fe = move( FEBase::build(dim, fe_type) );

  // Setup gauss quadrature
  qrule = QGauss( dim, fe_type.default_quadrature_order() );
  fe->attach_quadrature_rule (&qrule);

  // Enable calculations calculations
  fe->get_JxW();
  fe->get_phi();
  fe->get_dphi();

  // Jacobian
  for ( uint i=0; i<3; i++ )
  {
    vector<DenseSubMatrix<Number>> kk;
    for ( uint j=0; j<3; j++ ) kk.emplace_back( Ke );
    Ke_var.push_back(kk);
  }

  // RHS Vector
  for ( uint i=0; i<3; i++ )
    Fe_var.emplace_back(Fe);
}

/**
 *  Init the DoF map and the element matrices.
 *  This function must be called before (or at the beginning of)
 *  the jacobian and residual.
 */
void ViscoPlasticMaterial::reinit( const Elem & elem )
{
  SCOPELOG(1);

  fe->reinit( &elem );

  const DofMap & dof_map = system.get_dof_map();

  dof_map.dof_indices (&elem, dof_indices);
  dof_map.dof_indices (&elem, dof_indices_var[0], 0);
  dof_map.dof_indices (&elem, dof_indices_var[1], 1);
  dof_map.dof_indices (&elem, dof_indices_var[2], 2);

  uint n_dofs = dof_indices.size();
  uint n_dofsv = dof_indices_var[0].size();

  Ke.resize (n_dofs, n_dofs);
  for (uint vi=0; vi<3; vi++)
  for (uint vj=0; vj<3; vj++)
    Ke_var[vi][vj].reposition ( vi*n_dofsv, vj*n_dofsv, n_dofsv, n_dofsv );

  Fe.resize (n_dofs);
  for (uint vi=0; vi<3; vi++)
    Fe_var[vi].reposition ( vi*n_dofsv, n_dofsv );

}

/**
 *     Builds the jacobian of the element and assembles in the global _jacobian_.
 */
void ViscoPlasticMaterial::jacobian (const NumericVector<Number> & soln, SparseMatrix<Number> & jacobian)
{
  SCOPELOG(1);
  const vector<Real> & JxW = fe->get_JxW();
  const vector<vector<Real>> & phi = fe->get_phi();
  const vector<vector<RealGradient>> & dphi = fe->get_dphi();


}
/**
 *     Builds the RHS of the element and assembles in the global _residual_.
 */
void ViscoPlasticMaterial::residual (const NumericVector<Number> & soln, NumericVector<Number> & residual)
{
  SCOPELOG(1);
  const vector<Real> & JxW = fe->get_JxW();
  const vector<vector<Real>> & phi = fe->get_phi();
  const vector<vector<RealGradient>> & dphi = fe->get_dphi();
}
