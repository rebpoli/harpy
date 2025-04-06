
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

////  This goes in the jacobian and residual functions
//
//  const vector<Real> & JxW = fe_u->get_JxW();
//  const vector<vector<Real>> & phi_u = fe_u->get_phi();
//  const vector<vector<RealGradient>> & dphi_u = fe_u->get_dphi();

