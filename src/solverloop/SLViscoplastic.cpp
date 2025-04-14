
#include "solverloop/SLViscoplastic.h"

#include "config/Config.h"

#include "libmesh/equation_systems.h"
#include "libmesh/mesh.h"


/**
 *
 * SolverLoop for viscoplastic workflow.
 *
 * This is a single mesh workflow. The viscoplastic solver is 
 * the master solver, to load the mesh and the equation systems.
 *
 * The thermal is the slave solver, to use the ES and mesh from 
 * the master.
 *
 */
SLViscoplastic::SLViscoplastic( const Timestep & ts_ ) :
  Solverloop(ts_), viscoplastic( "viscoplastic", ts ),
  thermal( viscoplastic, "thermal" ),
  stress( viscoplastic, "stress" )
{
  // Initialize the ES
  EquationSystems & es = viscoplastic.es;
  es.init();

  // Init the solvers (must be after the es.init)
  viscoplastic.init();
  thermal.init();
  stress.init();

  // Initialize thermal solver. 
  thermal.solve();

  // Updates the stress
  thermal.update_coupler( stress );
  viscoplastic.update_coupler( stress );

  // Updates VP 2x to make sure the initial old solution is also updated in the coupler.
  thermal.update_coupler( viscoplastic );
  thermal.update_coupler( viscoplastic );

  stress.update_coupler( viscoplastic );
}

/**
 *   Loop around the solvers (within a single timeste) 
 *   until convergence is reached.
 */
void SLViscoplastic::solve() 
{
  SCOPELOG(1);
  
  // Update the temperature, sync the coupler
  thermal.solve(); 

  // Run the viscoplastic
  thermal.update_coupler( viscoplastic );
  stress.update_coupler( viscoplastic );
  viscoplastic.solve();

  viscoplastic.update_coupler( stress );
  thermal.update_coupler( stress );
  stress.solve();

  viscoplastic.export_exo("viscoplastic");
}

