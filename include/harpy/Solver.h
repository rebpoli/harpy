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

    virtual void solve()
      { flog << "Must be defined in the child class."; }
    virtual void export_results( string basename )
      { flog << "Must be defined in the child class."; }

  protected:
    const Timestep & ts;

};
