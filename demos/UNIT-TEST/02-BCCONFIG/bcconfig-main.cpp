
#include "config/BCConfig.h"
#include "solver/BC.h"
#include "base/HarpyInit.h"
#include "harpy/MeshInit.h"

#include "libmesh/equation_systems.h"
#include "libmesh/mesh.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/transient_system.h"


/*
 *
 * THIS TESTCASE IS DEPRECATED!!!
 *
 */

int main (int argc, char ** argv)
{
  // Read configuration file
  HarpyInit init( argc, argv );

  BCConfig bcconf("poroelastic");
  ilog << bcconf;

  ilog << "Test reftime";
  for ( int t=0 ; t<120 ; t++ )
    ilog << "    " << t << " => " << bcconf.get_reftime( t ) ;

  /** Create a dummy system **/
  libMesh::Mesh mesh(init.comm());
  MeshInit _mi( mesh );

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
  bc.update( 100 );
  ilog << bc;

  return 0;
}
