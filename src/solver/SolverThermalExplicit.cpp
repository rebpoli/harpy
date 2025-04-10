
#include "solver/SolverThermalExplicit.h"
#include "config/ModelConfig.h"
#include "config/BCConfig.h"
#include "harpy/Timestep.h"
#include "harpy/Material.h"
#include "util/String.h"
#include "libmesh/mesh.h"
#include "util/OutputOperators.h"

using harpy_string::iequals;

/**
 *
 */
SolverThermalExplicit::SolverThermalExplicit( string name_, const Timestep & ts_ ) :
    Solver(name_, ts_), bc_config(MODEL->boundary_config)
{}

/**
 *   Builds a simple structuer from the boundary configuration relating the
 *   material to its current temperature.
 *
 *   This will be used in the update_coupler later.
 */
void SolverThermalExplicit::solve()
{
  temperature_by_material.clear();
  double reftime = bc_config.get_reftime( ts.time ) ;
  if ( ! bc_config.entry_by_time.count( reftime ) ) flog << "Inconsistency in BCConfig. All _reftime_ must have an entry in entry_by_time.";

  const BCConfig::TimeEntry & timeentry = bc_config.entry_by_time.at( reftime );
  for ( auto & [ material, dbcs ] : timeentry.domain_bcs )
  for ( auto & dbc : dbcs )
  {
    if ( ! iequals ( dbc.vname , "T" ) ) continue;
    temperature_by_material[ material ] = dbc.value;
  }
}

/**
 *   Initialize the coupler for the target solver.
 */
void SolverThermalExplicit::init_trg_coupler( Solver & trg_solver )
{
  SCOPELOG(1);

  // For now, the Initial condition is just t=-1
  // So we should be safe by running 2 dry updates so that the we have an identicalk
  // current and old solution for the first timestep
  solve();
  update_coupler( trg_solver ) ;  
  update_coupler( trg_solver ) ;  
}

/**
 *
 */
void SolverThermalExplicit::update_coupler( Solver & trg_solver )
{
  SCOPELOG(1);

  Coupler & coupler = trg_solver.coupler;

  MeshBase * mesh = trg_solver.get_mesh();

  for (const auto & elem : mesh->active_element_ptr_range()) 
  {
    if ( elem->processor_id() != mesh->processor_id() ) continue;  /// Only in the active processor. NOTE: if this requires an element search, that would be a collective task! The loops have to be synchronous.

    uint eid = elem->id();
                                                                   ///
    Material * mat = trg_solver.get_material( *elem );
    mat->reinit( *elem );
    string mname = mat->name;


    // Update temperature
    double temperature = 0;
    if ( ! temperature_by_material.count( mname ) ) wlog << "No temperature defined for material '" << mname << "' at time " << ts.time << " (timestep " << ts.t_step << "). Using 0.";
    temperature = temperature_by_material[ mname ];

    // Update beta_e
    double beta_e = 0; 
    if ( ! mat->config.beta_e ) wlog << "No beta_e defined for material '" << mname << "' at time " << ts.time << " (timestep " << ts.t_step << "). Using 0.";
    else beta_e = *(mat->config.beta_e);

    // Update beta_d
    double beta_d = 0; 
    if ( ! mat->config.beta_d ) wlog << "No beta_d defined for material '" << mname << "' at time " << ts.time << " (timestep " << ts.t_step << "). Using 0.";
    else beta_d = *(mat->config.beta_d);

    if ( ! coupler.count(eid) ) coupler.emplace(eid, ElemCoupler(eid));

    ElemCoupler & ec = coupler.at( eid );

    // Save old temperature and store the new one
    ec.dbl_params_old["T"] = ec.dbl_params["T"];
    ec.dbl_params["T"].clear();
    for ( uint qp=0 ; qp<mat->qrule.n_points() ; qp++ )
      ec.dbl_params["T"].push_back( temperature );

  }
}

/**
 *    DUMPERS
 */
ostream& operator<<(ostream& os, const SolverThermalExplicit & m)
{
  os << "SolverThermalExplicit:" << endl;
  os << "       Name:                      " << m.name << endl;
  os << "       temperature_by_material:   " << m.temperature_by_material << endl;
  os << "       beta_e_by_material:        " << m.beta_e_by_material << endl;
  os << "       beta_d_by_material:        " << m.beta_d_by_material << endl;
  return os;
}
