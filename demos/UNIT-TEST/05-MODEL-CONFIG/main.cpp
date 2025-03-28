
#include "config/ModelConfig.h"
#include "base/HarpyInit.h"

int main (int argc, char ** argv)
{
  // Read configuration file
  HarpyInit init( argc, argv );

  ilog << "******************************************************";
  ilog << "                READS A COMPLETE MODEL";
  ilog << "******************************************************";


  ModelConfig mc( "model" );
  ilog << mc;

  return 0;
}
