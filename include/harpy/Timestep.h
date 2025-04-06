#pragma once

#include "base/Global.h"
#include <queue>

class TimestepConfig;

/**
 *
 *
 */
class Timestep {
  public:
    Timestep();
    bool test_end() const;
    double next();

    TimestepConfig & config;
    uint t_step;
    double dt, time;

    queue<double> tsqueue;  

  private:

    friend ostream& operator<<(ostream& os, const Timestep & ts);
};

ostream& operator<<(ostream& os, const Timestep & ts);
