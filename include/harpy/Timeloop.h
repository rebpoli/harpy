#pragma once

#include "base/Global.h"
#include "base/HarpyInit.h"
#include "config/TimeloopConfig.h"

class Poroelastic;

class Timeloop {
  public:
    Timeloop();
    ~Timeloop();
    void main_loop();

  private:
    Timestep ts;
};

