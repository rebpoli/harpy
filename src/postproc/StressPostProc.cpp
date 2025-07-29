
#include "postproc/StressPostProc.h"

#include "config/MaterialConfig.h"
#include "util/OutputOperators.h"

#include "libmesh/string_to_enum.h"

#include "libmesh/system.h"
#include "libmesh/equation_systems.h"
#include "libmesh/explicit_system.h"

/**
 *
 */
StressPostProc::StressPostProc( Material * refmat , ExplicitSystem & sys_ ) :
  ExplicitMaterial( refmat, sys_ )
{
  SCOPELOG(1);

  name = sys_.name();

//  dlog(1) << config;
}

/**
 *
 */
StressPostProc::~StressPostProc()
{ SCOPELOG(1); }


/**
 *
 */
void StressPostProc::init_fem()
{
  SCOPELOG(1);
  uint vid = system.variable_number( "sigtotXX" );
  DofMap & dof_map = system.get_dof_map();
  FEType fe_type = dof_map.variable_type(vid);

  fe = move( FEBase::build(3, fe_type) );

  // Qrule is imported from the referncematerial
  qrule = refmat->qrule;
  fe->attach_quadrature_rule (&qrule);

  // Enable calculations calculations
  fe->get_JxW();
  fe->get_phi();
  fe->get_dphi();
  fe->get_xyz();
}
