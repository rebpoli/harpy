#pragma once

#include "base/Global.h"
#include "harpy/Solver.h"
#include "harpy/Material.h"
#include "solver/BC.h"

#include "libmesh/mesh.h"
#include "libmesh/equation_systems.h"
#include "libmesh/transient_system.h"
#include "libmesh/nonlinear_implicit_system.h"

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

class SolverViscoplasticTrial : public Solver, 
                                public NonlinearImplicitSystem::ComputeResidual,
                                public NonlinearImplicitSystem::ComputeJacobian

{
  public:
    SolverViscoplasticTrial( string name, const Timestep & ts_ );
    ~SolverViscoplasticTrial();

    void init_materials();
    Material * get_material( const Elem & elem, bool reinit=0 );

    void solve();

    // Interface to set up the trial system
    virtual void jacobian (const NumericVector<Number> & soln,
                           SparseMatrix<Number> & jacobian,
                           NonlinearImplicitSystem & /*sys*/);
    virtual void residual (const NumericVector<Number> & soln,
                           NumericVector<Number> & residual,
                           NonlinearImplicitSystem & /*sys*/);


  private:
    
    void load_mesh();
    void set_dirichlet_bcs();
    void set_scalar_bcs() ;
    void set_unassigned_scalars();

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
