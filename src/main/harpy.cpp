#include <fstream>
#include <string>

#include "base/Global.h"
#include "base/HarpyInit.h"
#include "config/Config.h"
#include "harpy/Timeloop.h"
#include "util/Stopwatch.h"

int main (int argc, char ** argv)
{
  Stopwatch SW("FULL RUN");
  HarpyInit init( argc, argv );
  Timeloop timeloop;
  return 0;
}
