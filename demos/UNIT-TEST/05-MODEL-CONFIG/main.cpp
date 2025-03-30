
#include "config/ModelConfig.h"
#include "base/HarpyInit.h"

int main (int argc, char ** argv)
{
  // Read configuration file
  HarpyInit init( argc, argv );

  ilog << "******************************************************";
  ilog << "                READS A COMPLETE MODEL";
  ilog << "******************************************************";

  MODEL = new ModelConfig( "model/" );
  ilog << *MODEL;

  return 0;
}
