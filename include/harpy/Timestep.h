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

    // Signalize that the next step must be cut for the next n timesteps
    void cut() { cut_next = 5; } 

    TimestepConfig & config;
    uint t_step;
    double dt, time;

    int cut_next;

    queue<double> tsqueue;  

  private:

    friend ostream& operator<<(ostream& os, const Timestep & ts);
};

ostream& operator<<(ostream& os, const Timestep & ts);
