
#include "config/ModelConfig.h"
#include "config/BCConfig.h"
#include "solver/viscoplastic/VPBC.h"
#include "harpy/HarpyInit.h"

#include "libmesh/equation_systems.h"
#include "libmesh/mesh.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/transient_system.h"

using namespace harpy;
using namespace config;
using namespace solver::viscoplastic;

/*
 *
 */

int main (int argc, char ** argv)
{
  // Read configuration file
  HarpyInit init( argc, argv );

  MODEL = new ModelConfig( "model/" );
  BCConfig & bcconf = MODEL->boundary_config;;
  ilog << bcconf;

  ilog << "Test reftime";
  for ( int t=0 ; t<30 ; t++ )
    ilog << "    " << t/10. << " => " << bcconf.get_reftime( t/10. ) ;

  /** Create a dummy system **/
  libMesh::Mesh mesh(init.comm());

  EquationSystems es(mesh);
  TransientLinearImplicitSystem & sys = es.add_system<TransientLinearImplicitSystem> ("poroelastic");
  sys.add_variable("UX", FIRST, LAGRANGE);
  sys.add_variable("UY", FIRST, LAGRANGE);
  sys.add_variable("UZ", FIRST, LAGRANGE);
  sys.add_variable("P", FIRST, LAGRANGE);
  sys.add_variable("tri_dx", FIRST, LAGRANGE);
  sys.add_variable("tri_dy", FIRST, LAGRANGE);
  sys.add_variable("tri_dz", FIRST, LAGRANGE);
  sys.add_variable("K0", FIRST, LAGRANGE);
  BC bc(sys);
  bc.update( 1.5 );
  ilog << bc;

  return 0;
}
