
#pragma once

#include "base/Global.h"

/**
 *
 *  This file offers utilities to manage the outpt file (run/).
 *
 */

class Timestep;

namespace harpy_dirmanager
{

  string exo_filename( string fname, const Timestep & ts );

}
