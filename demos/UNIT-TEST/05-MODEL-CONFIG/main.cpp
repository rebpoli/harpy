
#include "config/ModelConfig.h"
#include "base/HarpyInit.h"

int main (int argc, char ** argv)
{
  // Read configuration file
  HarpyInit init( argc, argv );

  ModelConfig mc("poroelastic");
  ilog << mc;

  return 0;
}
