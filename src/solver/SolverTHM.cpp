
#include "solver/SolverTHM.h"

#include "libmesh/elem.h"

#include "config/ModelConfig.h"
#include "base/HarpyInit.h"

/**
 *  Creates the system and the materials.
 */
SolverTHM::SolverTHM( string name_ ) : Solver(), name(name_),
                                       config (MODEL->solver_config( name ) ) , 
                                       mesh( *LIBMESH_COMMUNICATOR ),
                                       es( mesh )
{
  dlog(1) << "SolverConfig: " << *config;
}

/**
 *   Deletes the owned data structure
 */
SolverTHM::~SolverTHM() {
  for ( auto & [ sid, mat ] : material_by_sid ) delete( mat );
  material_by_sid.clear();
}

/**
 *
 */
void SolverTHM::init_materials()
{
  // Init materials
  map<string, MaterialConfig> & materials = MODEL->materials;
  dlog(1) << materials;

  // ensures creation of all materials to the current mesh (local elems only)
  MeshBase & mesh = es.get_mesh();
  for ( const auto & elem : mesh.active_local_element_ptr_range() )
    get_material( *elem );
}

/**
 *   Returns the material for a given element.
 *   Creates a material if not existing.
 */
Material * SolverTHM::get_material( const Elem & elem, bool reinit )
{
  UNUSED(reinit);   /// TODO remove
  uint sid = elem.subdomain_id();
  if  ( ! material_by_sid.count( sid ) ) 
  {
    MeshBase & mesh = es.get_mesh();
    string sname = mesh.subdomain_name( sid );
//    if ( ! config->
//    material_by_sid[sid] = Material::Factory( sid );
  }

  return 0;
}

/**
 *
 */
void SolverTHM::solve()
{

}
