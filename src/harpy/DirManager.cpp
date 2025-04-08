
#include "harpy/DirManager.h"

#include "util/String.h"
#include <iomanip>
#include "harpy/Timestep.h"

namespace harpy_dirmanager
{

/**
 *
 */
string exo_filename( string fname, const Timestep & ts )
{
  using namespace harpy_string;

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
