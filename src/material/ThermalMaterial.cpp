
#include "material/ThermalMaterial.h"
#include "harpy/Coupler.h"
#include "config/MaterialConfig.h"
#include "util/OutputOperators.h"

#include "libmesh/system.h"
#include "libmesh/equation_systems.h"
#include "libmesh/explicit_system.h"



/**
 *
 */
ThermalMaterial::ThermalMaterial( suint sid_,
                                  const MaterialConfig & config, 
                                  ExplicitSystem & sys_ ) :
  MaterialExplicit( sid_, config, sys_ )
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
