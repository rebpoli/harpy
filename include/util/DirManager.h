
#pragma once

#include "harpy/Global.h"

/**
 *
 *  This file offers utilities to manage the outpt file (run/).
 *
 */

namespace timeloop { class Timestep; }

namespace util {

  using timeloop::Timestep;
  string exo_filename( string fname, const Timestep & ts );

}
