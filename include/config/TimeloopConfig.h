#pragma once

#include "base/Global.h"
#include <queue>

namespace libMesh { class System; }
class Timestep;

/**
 *
 *
 */
class TimeloopConfig {
  public:
    TimeloopConfig();

    uint n_steps;
    double dt, dt_k, dt_max, t_max, dt_min;
    vector<double> force_ts;

  private:

    queue<double> tsqueue;
    void build_tsvec();

    friend Timestep;
};

