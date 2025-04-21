#include <fstream>
#include <string>

#include "base/Global.h"
#include "base/HarpyInit.h"
#include "config/ModelConfig.h"
#include "harpy/Timeloop.h"
#include "util/Stopwatch.h"

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
