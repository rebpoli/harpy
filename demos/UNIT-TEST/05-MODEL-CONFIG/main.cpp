
#include "config/ModelConfig.h"
#include "config/InoutConfig.h"
#include "postproc/Probe.h"
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

  InoutConfig inout;
  ilog << inout;

  ProbeCol probes;
  ilog << probes;

  return 0;
}
