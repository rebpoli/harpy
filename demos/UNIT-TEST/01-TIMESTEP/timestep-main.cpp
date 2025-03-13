
#include "config/TimeloopConfig.h"
#include "base/HarpyInit.h"

int main (int argc, char ** argv)
{
  // Read configuration file
  HarpyInit init( argc, argv );

  ilog << "Starting Unit test";
  Timestep ts;
  while ( true  )
  {
    ilog << "Timestep " << ts.time << "s (dt=" << ts.dt << "s)";
    ts.next();
    if ( ts.test_end() ) break;
  }

  return 0;
}
