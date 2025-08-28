
#include "config/ModelConfig.h"
#include "config/InoutConfig.h"
#include "postproc/probes/Probe.h"
#include "harpy/HarpyInit.h"

using namespace harpy;
using namespace config;
using postproc::probes::ProbeCol;

int main (int argc, char ** argv)
{
  // Read configuration file
  HarpyInit init( argc, argv );

  ilog << "******************************************************";
  ilog << "                READS A COMPLETE MODEL";
  ilog << "******************************************************";

  MODEL = new ModelConfig( "model/" );
  ilog << *MODEL;

  InoutConfig inout;
  ilog << inout;

  ProbeCol probes;
  ilog << probes;

  return 0;
}
