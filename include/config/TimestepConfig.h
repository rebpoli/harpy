
#pragma once

#include "base/Global.h"
#include "util/String.h"

struct TimestepConfig
{
  TimestepConfig() :
    dt0(0.1), dtk(1), dt_max(1), t_max(100), max_steps(10), dt_min(1e-5) {};
  double dt0, dtk, dt_max, t_max, max_steps, dt_min;

  void set( string n, double v )
  {
    harpy_string::CIMap<double * > REFS = 
              { 
                { "dt0", &dt0 } ,
                { "dtk", &dtk } ,
                { "dt_min", &dt_min } ,
                { "tmax", &t_max } ,
                { "dtmax", &dt_max } ,
                { "max_steps", &max_steps } ,
              };

    if (! REFS.count( n ) ) flog << "Unknow Timestep key '" << n << "'. Check MODEL.";
    *(REFS[n]) = v;
  };
};

