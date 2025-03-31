
#include "solver/SolverTHM.h"

#include "libmesh/equation_systems.h"
#include "libmesh/elem.h"

#include "config/ModelConfig.h"

/**
 *  Creates the system and the materials.
 */
SolverTHM::SolverTHM( string name_ ) : Solver(), name(name_), es(0)
{
  // Read mesh, create EquationSystems
  SolverConfig * config = MODEL->solver_config( name );
  dlog(1) << "INIT SOLVER --- " << *config;
}

/**
 *   Deletes the owned data structure
 */
SolverTHM::~SolverTHM() { delete( es ) ; es = 0; }

/**
 *
 */
void SolverTHM::init_materials()
{
  // Init materials
  map<string, MaterialConfig> & materials = MODEL->materials;
  dlog(1) << materials;

  // ensures creation of all materials to the current mesh (local elems only)
  MeshBase & mesh = es->get_mesh();
  for ( const auto & elem : mesh.active_local_element_ptr_range() )
    get_material( *elem );
}

/**
 *   Returns the material for a given element.
 *   Creates a material if not existing.
 */
Material * SolverTHM::get_material( const Elem & elem )
{
  return 0;
}

/**
 *
 */
void SolverTHM::solve()
{

}
