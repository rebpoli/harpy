#pragma once

#include "base/Global.h"
#include <queue>

class Timestep;
namespace libMesh { class System; }

class TSCallback 
{
  public:
    virtual ~TSCallback() { }
    virtual void update( const Timestep & ts ) { 
      UNUSED(ts);
      flog<< "Should never be here."; 
    };
};

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

/**
 *
 *
 */
class Timestep {
  private:
    TimeloopConfig cfg;
    vector< TSCallback * > callbacks;
    void update() { for ( auto cb : callbacks ) cb->update( *this ); }

  public:
    Timestep() : cfg(), t_step(0), dt(cfg.dt), time(0) {}
    void add_callback( TSCallback & cb ) { callbacks.push_back( & cb ); cb.update( *this ); } 
    bool test_end() const;
    double next();
    uint t_step;
    double dt, time;

    friend ostream& operator<<(ostream& os, const Timestep & ts);
};

ostream& operator<<(ostream& os, const Timestep & ts);
