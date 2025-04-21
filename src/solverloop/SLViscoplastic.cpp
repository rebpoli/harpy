
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
  Solverloop(ts_), viscoplastic( "viscoplastic", ts ) , 
  thermal( viscoplastic, "thermal" )
{
  SCOPELOG(1);

  // Initialize the ES
  EquationSystems & es = viscoplastic.es;

  dlog(1) << "ES.init...";
  es.init();

  // Init the solvers (must be after the es.init)
  viscoplastic.init();
//  thermal.init();
}

/**
 *   Loop around the solvers (within a single timeste) 
 *   until convergence is reached.
 */
void SLViscoplastic::solve() 
{
  SCOPELOG(1);
  
  // Update the temperature, sync the coupler
//  thermal.solve(); 

  // Run the viscoplastic
  viscoplastic.solve();

  viscoplastic.export_exo("viscoplastic");
}

