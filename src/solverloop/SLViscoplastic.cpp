
#include "solverloop/SLViscoplastic.h"

#include "config/Config.h"

#include "libmesh/equation_systems.h"
#include "libmesh/mesh.h"


SLViscoplastic::SLViscoplastic( const Timestep & ts_ ) :
  Solverloop(ts_), viscoplastic( "viscoplastic", ts ), thermal( "thermal", ts )
{

  // Initialize structures
  thermal.init_coupler( viscoplastic );
}

/**
 *
 *
 */
void SLViscoplastic::solve() 
{
  SCOPELOG(1);
  
  // Update the temperature, sync the coupler
  thermal.solve(); 
  thermal.update_coupler( viscoplastic );

  // Run the viscoplastic
  viscoplastic.solve();
}

/**
 *
 *
 */
void SLViscoplastic::export_results() 
{
}
