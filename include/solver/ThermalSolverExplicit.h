#pragma once

#include "base/Global.h"
#include "harpy/Solver.h"
#include "util/String.h"
#include <map>

#include "postproc/ThermalPostProc.h"
#include "libmesh/explicit_system.h"
class BCConfig;

class ViscoplasticSolver;

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

class ThermalSolverExplicit : public Solver 
{
public:
  ThermalSolverExplicit( ViscoplasticSolver & ref_solver_, string name_ );
  virtual ~ThermalSolverExplicit() {};

  /// Init explicit materials (the ones to do the projections)
  void init_materials();

  /// The solver workers
  virtual void solve();

  /// A specialized system object
  ExplicitSystem & system;

private:
  ThermalPostProc * get_postproc( const Elem & elem );

  void project_to_system();
  void update_reference_solver();
  void setup_variables();

  /// Data organized from configuration
  CIMap<double> temperature_by_material;
  CIMap<double> initial_temperature_by_material;

  friend ostream& operator<<(ostream& os, const ThermalSolverExplicit & m);

  ViscoplasticSolver * ref_solver;

};

ostream& operator<<(ostream& os, const ThermalSolverExplicit & m);
