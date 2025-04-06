#pragma once

#include "base/Global.h"

class Timestep;

/**
 *
 *
 */

class Solver
{
  public:
    Solver( const Timestep & ts_ ) : ts(ts_) {};
    void solve();
    void export_results();

  protected:
    const Timestep & ts;
};
