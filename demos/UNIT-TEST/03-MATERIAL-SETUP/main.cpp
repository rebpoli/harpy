
#include "base/HarpyInit.h"
#include "config/ModelConfig.h"
#include "solver/ViscoplasticSolver.h"

#include "libmesh/mesh.h"
#include "libmesh/equation_systems.h"

#include "harpy/Timestep.h"

int main (int argc, char ** argv)
{
  // Read configuration file
  HarpyInit init( argc, argv );

  MODEL = new ModelConfig( "model/" );

  Timestep ts;
  ViscoplasticSolver st( "Viscoplastic", ts );

  return 0;
}
