
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
StressPostProc::StressPostProc( Material & refmat , ExplicitSystem & sys_ ) :
  ExplicitMaterial( refmat, sys_ )
{
  SCOPELOG(1);

  name = sys_.name();

//  dlog(1) << config;
  setup_variables();
}

/**
 *
 */
StressPostProc::~StressPostProc()
{ SCOPELOG(1); }


/**
 *   _sid_ is the subdomain id
 */
void StressPostProc::setup_variables()
{
  SCOPELOG(1);

  // Stress variable order is one under U
  if (! config.fem_by_var.count( "U" ) ) flog << "Undefined var setup for variable 'U'. Please revise model material.FEM section.";
  auto & femspec = config.fem_by_var.at("U");

  // TODO: fetch this info from configuration
  Order order = Utility::string_to_enum<Order>( femspec.order ) - 1;
  FEFamily fef = L2_LAGRANGE;
  if ( ! order ) fef = MONOMIAL;  // a constant is a monomial

  vector<string> sname = { "sigeff", "sigtot", "deviatoric", "plastic_strain", "plastic_strain_rate" };
  vector<string> sdir  = { "XX",  "YY",  "ZZ",  "XY",  "XZ",   "YZ" };

  set<subdomain_id_type> sids = { sid };
  for ( auto sn : sname ) 
  for ( auto sd : sdir )
    system.add_variable(sn+sd, order, fef, &sids);

  system.add_variable("von_mises", order, fef, &sids);
  system.add_variable("epskk", order, fef, &sids);
}

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
