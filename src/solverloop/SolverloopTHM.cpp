
#include "solverloop/SolverloopTHM.h"

#include "config/Config.h"

SolverloopTHM::SolverloopTHM( MeshBase & mesh, const Timestep & ts_ ) :
  es(mesh), ts(ts_), solver(es)
{
}

/**
 *
 *
 */
void SolverloopTHM::solve() 
{
  solver.solve();
}

/**
 *
 *
 */
void SolverloopTHM::export_results() 
{
}
