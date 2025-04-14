#pragma once

#include "base/Global.h"
#include "harpy/Solver.h"
#include "harpy/Material.h"
#include "solver/BC.h"

#include "libmesh/mesh.h"
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

    virtual void solve();

    // Interface to set up the trial system
    virtual void jacobian (const NumericVector<Number> & soln,
                           SparseMatrix<Number> & jacobian,
                           NonlinearImplicitSystem & /*sys*/);
    virtual void residual (const NumericVector<Number> & soln,
                           NumericVector<Number> & residual,
                           NonlinearImplicitSystem & /*sys*/);

    /// Updates the trg_solver coupler with info from this system
    void update_coupler( Solver & trg_solver );

    /// Updates the plastic strain in the coupler from info from the coupler itself
    void update_plastic_strain();

    /// A specialized system object
    TransientNonlinearImplicitSystem & system;

  private:
    
    void load_mesh();
    void set_dirichlet_bcs();
    void set_scalar_bcs() ;
    void set_unassigned_scalars();
    void init_materials();
    void add_scalar_vars();

    BC curr_bc;
};
