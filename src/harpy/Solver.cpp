#include "harpy/Solver.h"

#include "base/HarpyInit.h" // LIBMESH_COMMUNICATOR
#include "config/ModelConfig.h"
#include "harpy/Material.h"
#include "harpy/Timestep.h"
#include "harpy/DirManager.h"

#include "libmesh/system.h"
#include "libmesh/mesh.h"
#include "libmesh/elem.h"

/**
 *
 */
Solver::Solver( string name_, const Timestep & ts_ ) :
  name(name_), ts(ts_), own_es(1),
  config (MODEL->solver_config( name ) ),
  bc_config ( MODEL->boundary_config ),
  es( *(new EquationSystems( *(new Mesh(*LIBMESH_COMMUNICATOR)) ) ) ),
  coupler()
{}

/**
 *   When we want to share the ES with another solver (same mesh)
 */
Solver::Solver( EquationSystems & es_, string name_, const Timestep & ts_ ) :
  name(name_), ts(ts_), own_es(0),
  config (MODEL->solver_config( name ) ),
  bc_config ( MODEL->boundary_config ),
  es( es_ ),
  coupler()
{}

/**
 *   Deletes the equationsystem if it is owned by this object.
 *   Deletes the owned data structure.
 */
Solver::~Solver() {
  if ( own_es ) delete( &es ); 

  for ( auto & [ sid, mat ] : material_by_sid ) delete( mat );
  material_by_sid.clear();
}


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
    string sname = get_mesh().subdomain_name( sid );
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
  MeshBase & mesh = get_mesh();
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

/**
 *
 *
 */
void Solver::export_exo( string fn )
{
  using namespace harpy_dirmanager;
  fn = exo_filename( fn, ts );

  ExodusII_IO exo(get_mesh());
  exo.write_timestep_discontinuous ( fn, es, 1, ts.time );
}
