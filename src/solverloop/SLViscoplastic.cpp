
#include "solverloop/SLViscoplastic.h"

#include "config/Config.h"

#include "libmesh/equation_systems.h"
#include "libmesh/mesh.h"


SLViscoplastic::SLViscoplastic( const Timestep & ts_ ) :
  Solverloop(ts_), viscoplastic( "viscoplastic" )
{
  // Create solver and mesh 
}

/**
 *
 *
 */
void SLViscoplastic::solve() 
{
  viscoplastic.solve();
}

/**
 *
 *
 */
void SLViscoplastic::export_results() 
{
}
