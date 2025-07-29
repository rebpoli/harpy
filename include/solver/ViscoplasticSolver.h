#pragma once

#include "base/Global.h"

#include "config/InoutConfig.h"

#include "harpy/Solver.h"
#include "harpy/Material.h"
#include "postproc/ViscoplasticReport.h"
#include "material/ViscoPlasticMaterial.h"
#include "solver/BC.h"

#include "base/HarpyInit.h"

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
                                public NonlinearImplicitSystem::ComputeResidualandJacobian,
                                public System::Constraint
{
  public:
    ViscoplasticSolver( string name, Timestep & ts_ );
    virtual ~ViscoplasticSolver();

    /// Solution workflow
    virtual void solve();
    bool update_adaptive_timestep();
    void do_ts_cut();

    void constrain();

    /// Builds the linearized system of equations
    virtual void residual_and_jacobian (const NumericVector<Number> & soln,
                                        NumericVector<Number> * residual,
                                        SparseMatrix<Number> * jacobian, NonlinearImplicitSystem  & sys);


    /// Updates the plastic strain in the coupler from info from the coupler itself
    void update_plastic_strain();

    /// Helper
    ViscoPlasticMaterial * get_vp_material( const Elem & elem ) { return dynamic_cast<ViscoPlasticMaterial *> ( get_material(elem) ); }
    ViscoPlasticMaterial * get_vp_material( uint sid ) { return dynamic_cast<ViscoPlasticMaterial *> ( get_material(sid) ); }

    /// A specialized system object
    TransientNonlinearImplicitSystem & system;

    /// A system to output intermediate and posproc variables
    ExplicitSystem & stress_system;

  private:

    /// Material properties loaded once from configuration
    void load_mesh();
    void setup_variables();

    void set_dirichlet_bcs();
    void set_scalar_bcs() ;
    void set_unassigned_scalars();
    void init_materials();
    void add_scalar_vars();

    // Post processing
    void posproc_stresses();
    void export_results();
    ViscoplasticReport report;

    /// Boundary conditions for the current timetep
    BC curr_bc;

    unique_ptr<NumericVector<Number>> old_sol;
};
