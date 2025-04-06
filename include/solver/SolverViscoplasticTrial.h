#pragma once

#include "base/Global.h"
#include "harpy/Solver.h"
#include "harpy/Material.h"
#include "solver/BC.h"

#include "libmesh/mesh.h"
#include "libmesh/equation_systems.h"
#include "libmesh/transient_system.h"

#include <map>

class SolverConfig;
class Timestep;

namespace libMesh {
  class Elem;
  class MeshBase; 
}

/**
 *
 * Solution for NON - linear THM systems.
 * The underlying materials control which variables are
 * added into the system.
 *
 * The solver name controls the set of configurations
 * that are input.
 *
 * Ex: solver_name = 'poroelastic'
 *                 = 'thermal'
 *
 * This is always nonlinear. The materials
 * must be able to supply a jacobian and a residual for the
 * assigned element.
 *
 */

using namespace libMesh;

class SolverViscoplasticTrial : public Solver
{
  public:
    SolverViscoplasticTrial( string name, const Timestep & ts_ );
    ~SolverViscoplasticTrial();

    void init_materials();
    Material * get_material( const Elem & elem, bool reinit=0 );

    void solve();


  private:
    
    void load_mesh();
    void set_dirichlet_bcs();
    void add_scalar_vars();

    string name;

    // In order of initialization
    SolverConfig * config;
    Mesh mesh;
    EquationSystems es;

    TransientNonlinearImplicitSystem & system;
    BC curr_bc;

    map< uint, Material * > material_by_sid;
};
