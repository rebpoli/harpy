#pragma once

#include "harpy/Global.h"
#include "harpy/HarpyInit.h"
#include "timeloop/Timestep.h"

#include  "libmesh/mesh.h"

namespace timeloop {

using namespace libMesh;

class Timeloop {
  public:
    Timeloop();
    ~Timeloop();
    void main_loop();

  private:
    Timestep ts;
};

} //ns
