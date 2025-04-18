#include <fstream>
#include <string>

#include "base/Global.h"
#include "base/HarpyInit.h"
#include "config/ModelConfig.h"
#include "harpy/Timeloop.h"
#include "util/Stopwatch.h"

#include "libmesh/getpot.h"

void usage_error()
{
  ilog << "-------------------------------------------------";
  ilog << "Options: msh2exo ";
  ilog << " --in    filename  input mesh file";
  ilog << " --out   filename  output mesh file";
  ilog << "-------------------------------------------------";

  exit(1);
}

/**
 *
 */
template <typename T>
T assert_argument (GetPot & cl,
                   const std::string & argname,
                   const char * progname,
                   const T & defaultarg)
{
  if (!cl.search(argname))
    {
      elog << ("No " + argname + " argument found!") << std::endl;
      usage_error();
    }
  return cl.next(defaultarg);
}

/**
 *
 */
int main (int argc, char ** argv)
{
  Stopwatch SW("FULL RUN");
  HarpyInit init( argc, argv );

  /** **/
  GetPot cl(argc, argv);
  const std::string msh_fn = assert_argument(cl, "--in", argv[0], std::string(""));
  const std::string exo_fn = assert_argument(cl, "--out", argv[0], std::string(""));

  /** **/
  Mesh mesh(init.comm(), 3);
  ilog << "Reading mesh '" << msh_fn << "'...";
  mesh.read( msh_fn );
  mesh.all_second_order();

  /** **/
  ilog << "Exprting mesh '" << exo_fn << "'...";
  mesh.write( exo_fn );
  
  return 0;
}
