
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

  SolverTHM st( "Poroelastic" );

  return 0;
}
