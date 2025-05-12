#include "harpy/Solver.h"

#include "base/HarpyInit.h" // LIBMESH_COMMUNICATOR
#include "config/ModelConfig.h"
#include "harpy/Material.h"
#include "harpy/ExplicitMaterial.h"
#include "harpy/Timestep.h"
#include "harpy/DirManager.h"

#include "libmesh/system.h"
#include "libmesh/mesh.h"
#include "libmesh/elem.h"

/**
 *
 */
Solver::Solver( string name_, Timestep & ts_ ) :
  name(name_), ts(ts_), own_es(1),
  config (MODEL->solver_config( name_ ) ),
  bc_config ( MODEL->boundary_config ),
  mesh( *(new Mesh(*LIBMESH_COMMUNICATOR)) ),
  es( *(new EquationSystems(mesh) ) )
{}

/**
 *   When we want to share the ES with another solver (same mesh)
 */
Solver::Solver( Solver & ref, string name_ ) :
  name(name_), ts(ref.ts), own_es(0),
  config ( MODEL->solver_config( name_ ) ),
  bc_config ( MODEL->boundary_config ),
  mesh( ref.es.get_mesh() ),
  es( ref.es )
{
SCOPELOG(1);
}

/**
 *   Deletes the equationsystem if it is owned by this object.
 *   Deletes the owned data structure.
 */
Solver::~Solver() {
  SCOPELOG(1);

  if ( own_es ) 
  {
    delete( &es ); 
    delete( &mesh ) ;
  }

  for ( auto & [ sid, mat ] : material_by_sid ) delete( mat );
  material_by_sid.clear();
}

/*
 *
 */
void Solver::init()
{
  SCOPELOG(1);
  // Init FEM of the materials
  for ( auto & [ sid, mat ] : material_by_sid )
    mat->init_fem();
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
