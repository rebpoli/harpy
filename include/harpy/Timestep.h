#pragma once

#include "base/Global.h"
#include "config/TimeloopConfig.h"

/**
 *
 *
 */
class Timestep {
  private:
    TimeloopConfig cfg;

  public:
    Timestep() : cfg(), t_step(0), dt(cfg.dt), time(0) {}
    bool test_end() const;
    double next();
    uint t_step;
    double dt, time;

    friend ostream& operator<<(ostream& os, const Timestep & ts);
};

ostream& operator<<(ostream& os, const Timestep & ts);
