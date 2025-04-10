#include "harpy/Solver.h"
#include "config/ModelConfig.h"
#include "base/HarpyInit.h" // LIBMESH_COMMUNICATOR
#include "libmesh/elem.h"

/**
 *
 */
Solver::Solver( string name_, const Timestep & ts_ ) :
  name(name_), ts(ts_),
  config (MODEL->solver_config( name ) ),
  mesh( *LIBMESH_COMMUNICATOR ),
  es( mesh ),
  coupler() 
{}

/**
 *   Returns the material for a given element.
 *   Fails if not existing.
 */
Material * Solver::get_material( const Elem & elem )
{
  uint sid = elem.subdomain_id();
  // Consistency check
  if  ( ! material_by_sid.count( sid ) ) 
  {
    string sname = mesh.subdomain_name( sid );
    flog << "Cannot find material for SID '" << sname << "' (" << sid << ")";
  }

  Material * mat = material_by_sid.at(sid);

  return mat;
}

