#pragma once

#include "harpy/Global.h"
#include "solver/common/Solver.h"
#include "util/String.h"
#include <map>

#include "solver/common/ExplicitMaterial.h"
#include "libmesh/explicit_system.h"
#include "libmesh/point.h"

// FWD
namespace config { class BCConfig; }
namespace solver { namespace viscoplastic { class ViscoplasticSolver; } }
namespace util { class GridFile; }
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
namespace pressure {

using namespace libMesh;
using namespace solver::common;
using util::CIMap;
using util::GridFile;

/** **/
class PressureSolverFromFile : public Solver 
{
public:
  PressureSolverFromFile( ViscoplasticSolver * ref_solver_, string name_ );
  virtual ~PressureSolverFromFile();

  /// The solver workers
  virtual void solve();

private:

  ExplicitPressureMaterial * get_material( const Elem & elem );
  ExplicitPressureMaterial * get_material( uint sid );

  void init_materials();
  void init();

  void project_to_system();
  void update_reference_solver();
  void setup_variables();
  void read_from_file();

  friend ostream& operator<<(ostream& os, const PressureSolverFromFile & m);

  ViscoplasticSolver * ref_solver;
  ExplicitSystem & system;
  map< uint, ExplicitPressureMaterial * > material_by_sid;

  GridFile * grid;
  Point grid_origin;
};

ostream& operator<<(ostream& os, const PressureSolverFromFile & m);

}} // ns
