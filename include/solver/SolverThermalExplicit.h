#pragma once

#include "base/Global.h"
#include "harpy/Solver.h"
#include "util/String.h"
#include <map>
class BCConfig;

/**
 *
 *  This class provides coupling functions to a simple temperature 
 *  solver.
 *
 *  The temperature must be set explicitly in the model.
 *
 */

using harpy_string::CIMap;

class SolverThermalExplicit : public Solver 
{
public:
  SolverThermalExplicit( string name_, const Timestep & ts_ );

  virtual void solve();

  virtual void init_trg_coupler( Solver & trg_solver );
  virtual void update_coupler( Solver & trg_solver );

private:
  string name;
  BCConfig & bc_config;

  CIMap<double > temperature_by_material;
  CIMap<double > beta_e_by_material, beta_d_by_material;

  friend ostream& operator<<(ostream& os, const SolverThermalExplicit & m);
};

ostream& operator<<(ostream& os, const SolverThermalExplicit & m);
