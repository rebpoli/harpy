
#include "timeloop/Timestep.h"
#include "harpy/HarpyInit.h"
#include "config/ModelConfig.h"

using namespace harpy;
using namespace config;
using namespace timeloop;

int main (int argc, char ** argv)
{
  // Read configuration file
  HarpyInit init( argc, argv );

  MODEL = new ModelConfig( "model/" );
  
  dlog(1) << *MODEL;

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
