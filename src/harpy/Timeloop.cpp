#include "harpy/Timeloop.h"
#include "solverloop/SolverloopTHM.h"

#include "libmesh/equation_systems.h"

using namespace libMesh;

/**
 *
 * Starts the timeloop for the specified mesh.
 *
 */
Timeloop::Timeloop( MeshBase & mesh_ ) : mesh(mesh_), ts() 
{
  main_loop(); 
}

Timeloop::~Timeloop() { }

/**
 *
 * Main time loop. Controls the timestepping and callbacks associated with
 * the time. This is a control machine, not a worker.
 *
 */
void Timeloop::main_loop() 
{
  SolverloopTHM sloop( mesh, ts ); 

  while ( true )
  {
    ilog1 << "====================================================";
    ilog1 << "Solving timestep "<< ts.t_step<<" @ " << ts.time << "...";

    sloop.solve();
    sloop.export_results();

    // Avança os parametros do tstep
    ts.next();   
    if ( ts.test_end() ) break;
  }
}
