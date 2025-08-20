#include "harpy/Timeloop.h"
#include "restart/File.h"
#include "solverloop/SLViscoplastic.h"

using namespace libMesh;

/**
 *
 * Starts the timeloop.
 *
 */
Timeloop::Timeloop() : ts() 
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
  SLViscoplastic sloop( ts ); 

  // Read initial stress?
  if ( MODEL->sig0_file ) sloop.load_sig0_file( *(MODEL->sig0_file) );

  // TODO: Move as object member
  restart::File restart( MODEL->model_dir + "/restart.bin" );

  while ( true )
  {
    ilog1 << "====================================================";
    ilog1 << "Solving timestep "<< ts.t_step<<" @ " << ts.time << "s (dt=" << ts.dt << ") ...";

    sloop.solve();

    // TODO: we are writing a new restart at every timestep. Need to optimize.
    restart.write( sloop );

    ts.next();   
    if ( ts.test_end() ) break;
  }
}
