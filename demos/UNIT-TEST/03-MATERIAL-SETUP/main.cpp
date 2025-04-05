
#include "base/HarpyInit.h"
#include "config/ModelConfig.h"
#include "solver/SolverViscoplasticTrial.h"

#include "libmesh/mesh.h"
#include "libmesh/equation_systems.h"

#

int main (int argc, char ** argv)
{
  // Read configuration file
  HarpyInit init( argc, argv );

  MODEL = new ModelConfig( "model/" );

  SolverViscoplasticTrial st( "Viscoplastic" );

  return 0;
}
