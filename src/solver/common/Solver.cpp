#include "solver/common/Solver.h"

#include "harpy/HarpyInit.h" // LIBMESH_COMMUNICATOR
#include "config/ModelConfig.h"
#include "solver/common/ExplicitMaterial.h"
#include "timeloop/Timestep.h"
#include "util/DirManager.h"

#include "libmesh/system.h"
#include "libmesh/mesh.h"
#include "libmesh/elem.h"

//// Import children
#include "solver/thermal/TSolverCnt.h"
#include "solver/thermal/TSolverFile.h"
#include "solver/pressure/PSolverFile.h"

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

///**
// *   Returns the material for a given element.
// *   Fails if not existing.
// */
//Material * Solver::get_material( const Elem & elem )
//{
//  uint sid = elem.subdomain_id();
//  return get_material( sid );
//}
///**
// *   Returns a material by subdomain id
// */
//Material * Solver::get_material( uint sid )
//{
//  // Consistency check
//  if  ( ! material_by_sid.count( sid ) ) 
//  {
//    string sname = get_mesh().subdomain_name( sid );
//    flog << "Cannot find material for SID '" << sname << "' (" << sid << ")";
//  }

//  return material_by_sid.at(sid);
//}

/** 
 *
 * FACTORY
 *
 */

/** **/
Solver * SolverFactory::new_thermal( ViscoplasticSolver * ref )
{
  using namespace solver::thermal;

  string sname = "thermal";
  SolverConfig * config = MODEL->solver_config( sname );

  if ( config->external_file.is_defined() ) 
    return new ThermalSolverFromFile( ref,  sname );
  
  return new ThermalSolverConstant( ref,  sname );
}

/** **/
Solver * SolverFactory::new_pressure( ViscoplasticSolver * ref )
{
  using namespace solver::pressure;

  string sname = "pressure";
  SolverConfig * config = MODEL->solver_config( sname );

  if ( config->external_file.is_defined() ) 
    return new PressureSolverFromFile( ref,  sname );
  
  flog << "Unknown pressure solver!";
  return new PressureSolverFromFile( ref,  sname );  // Dummy line, to keep consistency
//  return new PressureSolverConstant( ref,  sname );
}

}} // ns
