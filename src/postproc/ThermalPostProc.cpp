
#include "postproc/ThermalPostProc.h"
#include "material/ViscoPlasticMaterial.h"
#include "config/MaterialConfig.h"
#include "util/OutputOperators.h"

#include "libmesh/system.h"
#include "libmesh/equation_systems.h"
#include "libmesh/explicit_system.h"



/**
 *
 */
ThermalPostProc::ThermalPostProc( ViscoPlasticMaterial & ref_material, ExplicitSystem & sys_ ) :
  ExplicitMaterial( ref_material, sys_ )
{ 
}


/**
 *
 */
void ThermalPostProc::init_fem()
{
  SCOPELOG(1);
  system.print_info();
  uint vid = system.variable_number( "T" );
  DofMap & dof_map = system.get_dof_map();
  FEType fe_type = dof_map.variable_type(vid);

  fe = move( FEBase::build(3, fe_type) );

  // Inherits the qrule from the reference material to
  // calculate in the same qp's
  if ( refmat )
    qrule = refmat->qrule;
  else
    qrule = QGauss( 3, fe_type.default_quadrature_order() );

  fe->attach_quadrature_rule (&qrule);

  // Enable calculations calculations
  fe->get_JxW();
  fe->get_phi();
  fe->get_dphi();
  fe->get_xyz();
}

