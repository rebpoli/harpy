#pragma once

#include "harpy/Global.h"
#include "solver/common/Solver.h"
#include "util/String.h"
#include <map>

#include "solver/common/ExplicitMaterial.h"
#include "libmesh/explicit_system.h"

namespace libMesh { class ExplicitSystem; }

namespace config { class BCConfig; }
namespace solver { namespace viscoplastic { class ViscoplasticSolver; } }
namespace util { class GridFile; }

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

using namespace util;
using namespace solver::common;

class ThermalSolverFromFile : public Solver 
{
public:
  ThermalSolverFromFile( ViscoplasticSolver * ref_solver_, string name_ );
  virtual ~ThermalSolverFromFile();

  /// The solver workers
  virtual void solve();


private:
  map< uint, ExplicitThermalMaterial * > material_by_sid;

  ExplicitThermalMaterial * get_material( const Elem & elem );
  ExplicitThermalMaterial * get_material( uint sid );

  /// Init explicit materials (the ones to do the projections)
  void init_materials();
  void init();

  void project_to_system();
  void update_reference_solver();
  void setup_variables();
  void read_from_file();


  friend ostream& operator<<(ostream& os, const ThermalSolverFromFile & m);

  ViscoplasticSolver * ref_solver;
  ExplicitSystem & system;

  GridFile * grid;
  libMesh::Point grid_origin;
};

ostream& operator<<(ostream& os, const ThermalSolverFromFile & m);

}} // ns
