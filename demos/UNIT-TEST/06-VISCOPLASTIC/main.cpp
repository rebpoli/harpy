
#include "config/ModelConfig.h"
#include "harpy/HarpyInit.h"
#include "timeloop/Timeloop.h"

using namespace harpy;
using namespace config;
using namespace timeloop;

int main (int argc, char ** argv)
{
  // Read configuration file
  HarpyInit init( argc, argv );

  ilog << "******************************************************";
  ilog << "                READS A COMPLETE MODEL";
  ilog << "******************************************************";

  MODEL = new ModelConfig( "model/" );
  ilog << *MODEL;

  // RUN!
  Timeloop main_loop;

  return 0;
}
