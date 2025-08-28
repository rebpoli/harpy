
#include "util/DirManager.h"

#include "util/String.h"
#include <iomanip>
#include "timeloop/Timestep.h"

namespace util
{

/**
 *
 */
string exo_filename( string fname, const Timestep & ts )
{
  using namespace util;

  string ret = "";

  replace_all(fname, "-", "_") ; 
  replace_all(fname, ".e", "") ; 

  fname = "run/exo/" + fname + string(".e");

  ostringstream ss;
  ss << setw(3) << setfill('0') << ts.t_step;
  fname = fname + string("-s.") + ss.str();

  return fname;
}


} // Namespace
