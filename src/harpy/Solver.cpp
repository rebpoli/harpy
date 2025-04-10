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
 *    This function needs the FEBase member to be initialized in the Material to fetch the xyz.
 */
void Solver::init_coupler()
{
  SCOPELOG(1);
  MeshBase & mesh = es.get_mesh();
  for ( const auto & elem : mesh.active_local_element_ptr_range() )
  {
    Material * mat = get_material( *elem );

    // Calc xyz
    const std::vector<Point> & xyz = mat->fe->get_xyz();
    mat->fe->reinit( elem );

    // Create element in the coupler
    uint eid = elem->id();
    if ( ! coupler.count(eid) ) coupler.emplace(eid, ElemCoupler(eid));
    ElemCoupler & ec = coupler.at( eid );

    // Feed the coupler
    for ( auto & pname : mat->required_material_properties )
    {
      const MaterialConfig & mconf = mat->config;
      mat->config.get_property( ec.dbl_params[pname], pname, xyz );
    }
  }
}
