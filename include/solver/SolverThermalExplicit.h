#pragma once

#include "base/Global.h"
#include "harpy/Solver.h"
#include "util/String.h"
#include <map>

#include "libmesh/explicit_system.h"
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
namespace libMesh { class ExplicitSystem; }

class SolverThermalExplicit : public Solver 
{
public:
  SolverThermalExplicit( Solver & ref, string name_ );

  virtual void solve();

  virtual void init_trg_coupler( Solver & trg_solver );
  virtual void update_coupler( Solver & trg_solver );
  void init_materials();

  /// A specialized system object
  ExplicitSystem & system;

private:

  CIMap<double> temperature_by_material;

  friend ostream& operator<<(ostream& os, const SolverThermalExplicit & m);

};

ostream& operator<<(ostream& os, const SolverThermalExplicit & m);
