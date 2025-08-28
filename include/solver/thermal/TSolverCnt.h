#pragma once

#include "harpy/Global.h"
#include "solver/common/Solver.h"
#include "util/String.h"
#include <map>
#include "solver/common/ExplicitMaterial.h"
#include "libmesh/explicit_system.h"

class BCConfig;
class ViscoplasticSolver;

namespace libMesh { class ExplicitSystem; }

/**
 *
 *  This class provides coupling functions to a simple temperature 
 *  solver.
 *
 *  The temperature must be set explicitly in the model to a constant
 *  at every timestep, for each material.
 *
 *  Each material has a single temperature value at each timestep.
 *
 */


namespace solver {
namespace thermal {

using namespace solver::common;
using util::CIMap;

class ThermalSolverConstant : public Solver 
{
public:
  ThermalSolverConstant( ViscoplasticSolver * ref_solver_, string name_ );
  virtual ~ThermalSolverConstant() ;

  /// The solver workers
  virtual void solve();

private:
  map< uint, ExplicitThermalMaterial * > material_by_sid;

  ExplicitThermalMaterial * get_material( const Elem & elem );
  ExplicitThermalMaterial * get_material( uint sid );

  void init();
  void init_materials();

  void project_to_system();
  void update_reference_solver();
  void setup_variables();

  /// Data organized from configuration
  map<uint, double> temperature_by_sid;
  map<uint, double> initial_temperature_by_sid;

  ViscoplasticSolver * ref_solver;
  ExplicitSystem & system;

  friend ostream& operator<<(ostream& os, const ThermalSolverConstant & m);
};

ostream& operator<<(ostream& os, const ThermalSolverConstant & m);

}} // ns
