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
class ViscoPlasticMaterial;

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

class ViscoplasticSolver : public Solver, 
                                public NonlinearImplicitSystem::ComputeResidualandJacobian
{
  public:
    ViscoplasticSolver( string name, const Timestep & ts_ );
    virtual ~ViscoplasticSolver();

    /// Solution workflow
    virtual void solve();

    /// Builds the linearized system of equations
    virtual void residual_and_jacobian (const NumericVector<Number> & soln,
                                        NumericVector<Number> * residual,
                                        SparseMatrix<Number> * jacobian, NonlinearImplicitSystem  & sys);


    /// Updates the plastic strain in the coupler from info from the coupler itself
    void update_plastic_strain();

    /// A specialized system object
    TransientNonlinearImplicitSystem & system;

    /// A system to output intermediate and posproc variables
    ExplicitSystem & stress_system;

    /// The material can be provided for the dependent solvers
    ViscoPlasticMaterial * get_material( const Elem & elem );

  private:

    /// Material properties loaded once from configuration
    void load_mesh();

    void set_dirichlet_bcs();
    void set_scalar_bcs() ;
    void set_unassigned_scalars();
    void init_materials();
    void add_scalar_vars();
    void posproc_stresses();

    /// Boundary conditions for the current timetep
    BC curr_bc;
};
