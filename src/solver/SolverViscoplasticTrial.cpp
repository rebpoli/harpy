
#include "solver/SolverViscoplasticTrial.h"

#include "libmesh/elem.h"
#include "libmesh/boundary_info.h"

#include "config/ModelConfig.h"
#include "base/HarpyInit.h"
#include "util/MeshUtils.h"

/**
 *  Creates the system and the materials.
 *
 *  This object owns the mesh and the rquation system.
 */
SolverViscoplasticTrial::SolverViscoplasticTrial( string name_ ) : Solver(), name(name_),
                                       config (MODEL->solver_config( name ) ) , 
                                       mesh( *LIBMESH_COMMUNICATOR ),
                                       es( mesh )
{
  dlog(1) << "SolverViscoplasticTrial: " << *config;
  load_mesh();
  init_materials();
}

/**
 *    Reads and prepares the mesh for usage.
 */
void SolverViscoplasticTrial::load_mesh()
{
  string fn = config->mesh_filename;
  dlog(1) << "Reading mesh '" << fn << "'...";
  mesh.read( fn );
  dump_mesh( mesh );
}

/**
 *   Deletes the owned data structure
 */
SolverViscoplasticTrial::~SolverViscoplasticTrial() 
{
  for ( auto & [ sid, mat ] : material_by_sid ) delete( mat );
  material_by_sid.clear();
}

/**
 *
 */
void SolverViscoplasticTrial::init_materials()
{
  SCOPELOG(1);
  // ensures creation of all materials to the current mesh (local elems only)
  MeshBase & mesh = es.get_mesh();
  for ( const auto & elem : mesh.active_local_element_ptr_range() )
    get_material( *elem );
}

/**
 *   Returns the material for a given element.
 *   Creates a material if not existing.
 */
Material * SolverViscoplasticTrial::get_material( const Elem & elem, bool reinit )
{
  uint sid = elem.subdomain_id();
  if  ( ! material_by_sid.count( sid ) ) 
    material_by_sid[sid] = Material::Factory( sid, es.get_mesh(), *config );

  Material * mat = material_by_sid.at(sid);

  if ( reinit ) mat->reinit();

  return mat;
}

/**
 *
 */
void SolverViscoplasticTrial::solve()
{

}

