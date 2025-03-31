
#include "base/HarpyInit.h"
#include "config/ModelConfig.h"
#include "harpy/MeshInit.h"
#include "solver/SolverTHM.h"

#include "libmesh/mesh.h"
#include "libmesh/equation_systems.h"

#

int main (int argc, char ** argv)
{
  // Read configuration file
  HarpyInit init( argc, argv );

  MODEL = new ModelConfig( "model/" );

  /** Create a dummy system **/
  libMesh::Mesh mesh(init.comm());
//  MeshInit mi( mesh );

  EquationSystems es( mesh );
  SolverTHM st( es );

  return 0;
}
