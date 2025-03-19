
#include "config/BCConfig.h"
#include "base/HarpyInit.h"

int main (int argc, char ** argv)
{
  // Read configuration file
  HarpyInit init( argc, argv );

  BCConfig bc("poroelastic");
  ilog << bc;

  return 0;
}
