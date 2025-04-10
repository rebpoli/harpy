#include "harpy/Solver.h"

#include "base/HarpyInit.h" // LIBMESH_COMMUNICATOR
#include "config/ModelConfig.h"
#include "harpy/Material.h"
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

/**
 *    Fetches information from the configuration and feeds the object coupler.
 *
 */
void Solver::init_coupler()
{
  SCOPELOG(1);
  MeshBase & mesh = es.get_mesh();
  for ( const auto & elem : mesh.active_local_element_ptr_range() )
  {
    Material * mat = get_material( *elem );


    // Create element in the coupler
    uint eid = elem->id();
    if ( ! coupler.count(eid) ) coupler.emplace(eid, ElemCoupler(eid));
    ElemCoupler & ec = coupler.at( eid );
    mat->init_coupler( elem, ec );
  }
}
