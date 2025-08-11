#pragma once

#include "base/Global.h"
#include "harpy/Solver.h"
#include "util/String.h"
#include <map>

#include "harpy/ExplicitMaterial.h"
#include "libmesh/explicit_system.h"

class BCConfig;
class ViscoplasticSolver;
class GridFile;

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

using harpy_string::CIMap;
namespace libMesh { class ExplicitSystem; }

class PressureSolverFromFile : public Solver 
{
public:
  PressureSolverFromFile( Solver * ref_solver_, string name_ );
  virtual ~PressureSolverFromFile();

  /// The solver workers
  virtual void solve();

  /// A specialized system object
  ExplicitSystem & system;

private:
  ExplicitMaterial * get_explicit_material( const Elem & elem ) { return dynamic_cast<ExplicitMaterial *> ( get_material( elem ) ); }

  /// Init explicit materials (the ones to do the projections)
  void init_materials();

  void project_to_system();
  void update_reference_solver();
  void setup_variables();
  void read_from_file();

  /// Data holder
  GridFile * grid;

  /// Helper stuff
  libMesh::Point grid_origin;

  friend ostream& operator<<(ostream& os, const PressureSolverFromFile & m);


};

ostream& operator<<(ostream& os, const PressureSolverFromFile & m);
