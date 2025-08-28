
#include "solverloop/SLViscoplastic.h"

#include "config/Config.h"

#include "libmesh/equation_systems.h"
#include "libmesh/mesh.h"

#include "restart/File.h"

namespace solverloop {

using solver::common::SolverFactory;

/**
 *
 * SolverLoop for viscoplastic workflow.
 *
 * This is a single mesh workflow. The viscoplastic solver is 
 * the master solver, to load the mesh and the equation systems.
 *
 * The thermal and pressure are slave solvers, to use the ES and mesh from 
 * the master.
 *
 */
SLViscoplastic::SLViscoplastic( Timestep & ts_ ) :
  Solverloop(ts_), viscoplastic(0), thermal(0) , pressure(0)
{
  SCOPELOG(1);

  viscoplastic = new ViscoplasticSolver( "viscoplastic", ts );
  thermal = SolverFactory::new_thermal( viscoplastic );
  pressure = SolverFactory::new_pressure( viscoplastic );

  // Initialize the ES
  EquationSystems & es = viscoplastic->es;

  dlog(1) << "ES.init...";
  es.init();

  // Init the solvers (must be after the es.init)
  viscoplastic->init();
  thermal->init();
  pressure->init();

  // From the configuration setup, initialize the strain
  thermal->solve(); 
  pressure->solve(); 
}

/**
 *
 */
SLViscoplastic::~SLViscoplastic()
{
  delete(viscoplastic);
  delete(thermal);
  delete(pressure);
}

/**
 *   Loop around the solvers (within a single timeste) 
 *   until convergence is reached.
 */
void SLViscoplastic::solve() 
{
  SCOPELOG(1);
  
  // Update the temperature and pressure, sync the coupler
  thermal->solve(); 
  pressure->solve(); 

  // Run the viscoplastic
  viscoplastic->solve();

  viscoplastic->export_exo("viscoplastic");
}

/**
 *
 *
 */
void SLViscoplastic::load_sig0_file( string filename )
{
  SCOPELOG(1);
  restart::File restart( filename ) ;
  restart.read( *this );
}

} // ns
