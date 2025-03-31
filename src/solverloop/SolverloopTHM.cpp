
#include "solverloop/SolverloopTHM.h"

#include "config/Config.h"

#include "libmesh/equation_systems.h"
#include "libmesh/mesh.h"


SolverloopTHM::SolverloopTHM( const Timestep & ts_ ) :
  Solverloop(ts_), poroelastic( "Poroelastic" )
{
  // Create solver and mesh 
}

/**
 *
 *
 */
void SolverloopTHM::solve() 
{
  poroelastic.solve();
}

/**
 *
 *
 */
void SolverloopTHM::export_results() 
{
}
