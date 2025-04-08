
#include "config/ModelConfig.h"
#include "base/HarpyInit.h"
#include "harpy/Timeloop.h"

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
