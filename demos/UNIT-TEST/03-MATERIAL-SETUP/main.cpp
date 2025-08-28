
#include "harpy/HarpyInit.h"
#include "config/ModelConfig.h"
#include "solver/viscoplastic/VPSolver.h"

#include "libmesh/mesh.h"
#include "libmesh/equation_systems.h"

#include "timeloop/Timestep.h"

using namespace harpy;
using namespace solver::viscoplastic;
using namespace timeloop;
using namespace config;

int main (int argc, char ** argv)
{
  // Read configuration file
  HarpyInit init( argc, argv );

  MODEL = new ModelConfig( "model/" );

  Timestep ts;
  ViscoplasticSolver st( "Viscoplastic", ts );

  return 0;
}
