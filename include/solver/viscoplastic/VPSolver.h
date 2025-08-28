#pragma once

#include "harpy/Global.h"

#include "config/InoutConfig.h"

#include "solver/common/Solver.h"
#include "postproc/report/VPReport.h"
#include "solver/viscoplastic/VPMaterial.h"
#include "solver/viscoplastic/VPMatBC.h"
#include "solver/viscoplastic/VPMatEG.h"
#include "solver/viscoplastic/VPBC.h"

#include "harpy/HarpyInit.h"

#include "libmesh/mesh.h"
#include "libmesh/transient_system.h"
#include "libmesh/nonlinear_implicit_system.h"

#include <map>

namespace restart { class File; } 

namespace libMesh { class Elem; class MeshBase; }

namespace timeloop { class Timestep; }
namespace confib { class SolverConfig; }

namespace solver {
namespace viscoplastic {

using postproc::report::ViscoplasticReport;
using timeloop::Timestep;
using solver::common::Solver;

class ViscoPlasticMaterial; 
class ViscoPlasticMaterialBC; 


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
    virtual void init();

    bool update_adaptive_timestep();
    void do_ts_cut();

    void constrain();

    /// Builds the linearized system of equations
    virtual void residual_and_jacobian (const NumericVector<Number> & soln,
                                        NumericVector<Number> * residual,
                                        SparseMatrix<Number> * jacobian, NonlinearImplicitSystem  & sys);


    /// Updates the plastic strain in the coupler from info from the coupler itself
    void update_plastic_strain();

    // Helpers
    bool is_eg() { return system.has_variable( "UegX" ); }
    
    /// A specialized system object
    TransientNonlinearImplicitSystem & system;

    /// A system to output intermediate and posproc variables
    ExplicitSystem & stress_system;

    ///
    ViscoPlasticMaterial * get_material( uint sid );
    ViscoPlasticMaterialBC * get_mat_bc( uint sid );
    ViscoPlasticMaterial * get_material( const Elem & elem );
    ViscoPlasticMaterialBC * get_mat_bc( const Elem & elem );

  private:
    map< uint, ViscoPlasticMaterial *> material_by_sid;
    map< uint, ViscoPlasticMaterialBC *> matbc_by_sid;

    // EG engine for the interface discontinuities
    VPMatEG * vpmat_eg;

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

    friend class restart::File;
    friend ViscoplasticReport;
};

}} // ns
