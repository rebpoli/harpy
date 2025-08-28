#include "solver/common/Solver.h"

#include "harpy/HarpyInit.h" // LIBMESH_COMMUNICATOR
#include "config/ModelConfig.h"
#include "solver/common/ExplicitMaterial.h"
#include "timeloop/Timestep.h"
#include "util/DirManager.h"

#include "libmesh/system.h"
#include "libmesh/mesh.h"
#include "libmesh/elem.h"

namespace solver {
namespace common {

using config::MODEL;

/**
 *
 */
Solver::Solver( string name_, Timestep & ts_ ) :
  name(name_), ts(ts_), own_es(1),
  config (MODEL->solver_config( name_ ) ),
  bc_config ( MODEL->boundary_config ),
  es( *(new EquationSystems( *(new Mesh(*LIBMESH_COMMUNICATOR)) ) ) )
{ }

/**
 *   When we want to share the ES with another solver (same mesh)
 */
Solver::Solver( Solver * ref, string name_ ) :
  name(name_), ts(ref->ts), own_es(0),
  config ( MODEL->solver_config( name_ ) ),
  bc_config ( MODEL->boundary_config ),
  es( ref->es )
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
    MeshBase & mesh = es.get_mesh();
    delete( &es ); 
    delete( & mesh ) ;
  }

}

/**
 *
 *
 */
void Solver::export_exo( string fn )
{
  using util::exo_filename;
  fn = exo_filename( fn, ts );

  ExodusII_IO exo(get_mesh());
  exo.write_timestep_discontinuous ( fn, es, 1, ts.time );
}

}} // ns
