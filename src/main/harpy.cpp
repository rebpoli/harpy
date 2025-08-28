#include <fstream>
#include <string>

#include "harpy/Global.h"
#include "harpy/HarpyInit.h"
#include "config/ModelConfig.h"
#include "timeloop/Timeloop.h"
#include "util/Stopwatch.h"

using namespace util;
using namespace harpy;
using namespace timeloop;
using namespace config;

int main (int argc, char ** argv)
{
  Stopwatch SW("FULL RUN");
  HarpyInit init( argc, argv );

  // Load the model config
  MODEL = new ModelConfig( "model/" );

  // Starts the timeloop
  Timeloop timeloop;

  return 0;
}
