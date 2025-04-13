#pragma once

#include "base/Global.h"

#include "harpy/Solver.h"

class SolverStress : public Solver
{
public:
  SolverStress( Solver & ref , string name_ );

  ExplicitSystem & system;

  virtual void solve();
  virtual void init_materials();

  void update_coupler( Solver & trg_solver );
};

