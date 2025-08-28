

#include "solver/viscoplastic/VPSolver.h"
#include "config/ModelConfig.h"
#include "timeloop/Timestep.h"
#include "util/DirManager.h"
#include "postproc/stress/StressPostProc.h"
#include "util/MeshUtils.h"
#include "util/Messages.h"
#include "util/OutputOperators.h"
#include "util/Stopwatch.h"

#include "libmesh/boundary_info.h"
#include "libmesh/elem.h"
#include "libmesh/dof_map.h"
#include "libmesh/transient_system.h"
#include "libmesh/nonlinear_implicit_system.h"
#include "libmesh/nonlinear_solver.h"

#include "libmesh/dirichlet_boundaries.h"
#include "libmesh/const_function.h"

#include "libmesh/string_to_enum.h"
#include "libmesh/enum_fe_family.h"
#include "libmesh/enum_order.h"

//#include <petscsnes.h>
#include "libmesh/petsc_nonlinear_solver.h"

namespace solver {
namespace viscoplastic {

using config::MODEL;
using config::SolverConfig;
using namespace util;

/**
 *  Creates the system and the materials.
 *
 *  This object owns the mesh and the rquation system.
 */
ViscoplasticSolver::ViscoplasticSolver( string name_, Timestep & ts_ ) : 
                   Solver( name_, ts_ ), 
                   system( es.add_system<TransientNonlinearImplicitSystem> ( name ) ),
                   stress_system(es.add_system<ExplicitSystem> ( name+"-stress" )),
                   report(*this), curr_bc( system )
{
  SCOPELOG(1);
  dlog(1) << "ViscoplasticSolver: " << *config;

  // Init equationsystems flow
  load_mesh();
  init_materials();
  setup_variables();
  report.init();

  system.nonlinear_solver->residual_and_jacobian_object = this;
  system.attach_constraint_object( *this );
}

/**
 *
 */
ViscoplasticSolver::~ViscoplasticSolver() {}

/**  **/
void ViscoplasticSolver::init()
{
  SCOPELOG(1);
  for ( auto & [ sid, mat ] : material_by_sid )
    mat->init_fem();
}

/**
 *    Adds all variables needed for the solver and its children (stresses, for example)
 *    
 *    Note: Variables are added for the whole mesh. If a material do not use them,
 *          Dirichlet BCs should be added to all nodes to prevent their solution.
 */
void ViscoplasticSolver::setup_variables()
{
  // Add displacements variables to the current subdomain ID
  if (! config->fem_by_var.count( "U" ) ) flog << "Undefined var setup for variable 'U'. Please revise model material.FEM section.";
  auto & femspec = config->fem_by_var.at("U");

  {
    set<string> KNOWN_FAMILY = { "LAGRANGE", "ENRICHED_GALERKIN" };
    if ( ! KNOWN_FAMILY.count( femspec.family ) ) flog << "FEM family not supported for this solver '" << femspec.family << "'.";

    Order order = Utility::string_to_enum<Order>( femspec.order ) ;

//    dlog(1) << "Setting up variable 'U' for ViscoplasticSolver ...";
//    dlog(1) << "     Order:" << order;
//    dlog(1) << "     FEFamily:" << fe_family;

    system.add_variable( "UX", order, LAGRANGE );
    system.add_variable( "UY", order, LAGRANGE );
    system.add_variable( "UZ", order, LAGRANGE );

    if ( femspec.family == "ENRICHED_GALERKIN")
    {
      system.add_variable( "UegX", CONSTANT, MONOMIAL );
      system.add_variable( "UegY", CONSTANT, MONOMIAL );
      system.add_variable( "UegZ", CONSTANT, MONOMIAL );
    }
  }

  // Stresses
  {
    Order order = Utility::string_to_enum<Order>( femspec.order ) - 1;
    FEFamily fef = L2_LAGRANGE;
    if ( ! order ) fef = MONOMIAL;  // a constant is a monomial


    /* Vectors */

    // Invariants
    vector<string> sdir_vec  = { "X",  "Y",  "Z" };

    // S1, S2, S3
    for ( uint i=0 ; i<3; i++ )
    for ( auto sd : sdir_vec )
    {
      string vn = "S"+to_string(i+1)+sd;
      stress_system.add_variable( vn, order, fef );
      dlog(1) << "Added variable '" << vn << "'";
    }

    // Sterz1, Sterz2, Sterz3 (Note: X,Y,Z of a vector must be added in sequence)
    for ( uint i=0 ; i<3; i++ )
    for ( auto sd : sdir_vec )
    {
      string vn = "Sterz"+to_string(i+1)+sd;
      stress_system.add_variable( vn, order, fef );
    }

    for ( uint i=0 ; i<3; i++ )
      stress_system.add_variable( "S"+to_string(i+1)+"_mag" , order, fef);

    for ( uint i=0 ; i<3; i++ )
      stress_system.add_variable( "Sterz"+to_string(i+1)+"_mag" , order, fef);

    /* TENSORS */

    vector<string> sname = { "sigeff", "sigtot", "deviatoric", "plastic_strain", "plastic_strain_rate", "initial_stress" };
    vector<string> sdir  = { "XX",  "YY",  "ZZ",  "XY",  "XZ",   "YZ" };

    dlog(1) << "Adding variables in stress post proc ...";
    for ( auto sn : sname ) 
    for ( auto sd : sdir )
      stress_system.add_variable( sn+sd, order, fef );

    /* Scalars */
    stress_system.add_variable( "von_mises", order, fef);
    stress_system.add_variable( "epskk", order, fef );
    stress_system.add_variable( "F", order, fef );

    stress_system.add_variable( "S_invarQ" , order, fef);
    stress_system.add_variable( "S_invarP" , order, fef);
    stress_system.add_variable( "S_invarQ_div_P" , order, fef);
    stress_system.add_variable( "Sterz_invarQ" , order, fef);
    stress_system.add_variable( "Sterz_invarP" , order, fef);
    stress_system.add_variable( "Sterz_invarQ_div_P" , order, fef);

  }

  // Scalars and penalties
  add_scalar_vars();
}

/**
 *
 *
 */
void ViscoplasticSolver::constrain()
{
  SCOPELOG(1);
  set_scalar_bcs();
}

/**
 *   Creates all the needed materials for the solution (one per subdomain ID).
 *   Initializes the material_by_sid structure.
 */
void ViscoplasticSolver::init_materials()
{
  SCOPELOG(1);
  MeshBase & mesh = get_mesh();

  set<MaterialConfig> & materials = MODEL->materials;

  // ensures creation of all materials to the current mesh (local elems only)
  for ( const auto & elem : mesh.active_element_ptr_range() )
  {
    uint sid = elem->subdomain_id();
    if  ( material_by_sid.count( sid ) ) continue;

    string sname = mesh.subdomain_name( sid );
    SolverConfig & svr_config = *( this->config );
    if ( ! svr_config.mat_config_by_name.count( sname ) ) flog << "Cannot find material configuration by name for subdomain '" << sname << "'. The model is inconsistent.";
    auto & mat_conf_id = svr_config.mat_config_by_name.at( sname );

    // Build material object
    string mat_name = mat_conf_id.name, mat_cfg = mat_conf_id.cfg;
    MaterialConfig mckey( mat_name, mat_cfg );
    auto it = materials.find( mckey );
    if ( it == materials.end() ) flog << "Cannot find material description for '" << sname << "'. The model is inconsistent.";

    // This object has all physical properties (por, perm, alpha, ...)
    const MaterialConfig & mat_conf = *it;

    dlog(1) << "Resoved material:" << mat_conf;
    ViscoPlasticMaterial * vpmat =
      new ViscoPlasticMaterial( sid, &mat_conf, system, *this );
    material_by_sid[sid] = vpmat;
  }
}

//  return dynamic_cast<ViscoPlasticMaterial *> ( material_by_sid.at(sid) );

/**
 *
 */
void ViscoplasticSolver::add_scalar_vars()
{
  SCOPELOG(1);
  for ( auto & sv : MODEL->boundary_config.scalars )
  {
    dlog(1) << "Adding scalar var " << sv.name << " ...";
    system.add_variable(sv.name, FIRST, SCALAR);
  }
}

/**
 *    Reads and prepares the mesh for usage.
 */
void ViscoplasticSolver::load_mesh()
{
  string fn = config->mesh_filename;
  dlog(1) << "Reading mesh '" << fn << "'...";
  MeshBase & mesh = get_mesh();
  mesh.read( fn );

  // Wee need a second order mesh 
  {
    Stopwatch sw("mesh.all_second_order()");
    mesh.all_second_order();
  }
}

/**
 *
 */
void ViscoplasticSolver::set_dirichlet_bcs()
{
  SCOPELOG1(1);

  const MeshBase & mesh = get_mesh();
  const BoundaryInfo & bi = mesh.get_boundary_info();
  DofMap & dof_map = system.get_dof_map();
  dof_map.get_dirichlet_boundaries()->clear();

  // The dirichlet structure in curr_bc is already processed and can be added easily
  for ( auto & dbc : curr_bc.dirichlet ) 
  {
    dlog(1) << "["<< fmt_i(ts.t_step) << "] Adding dirichlet boundary condition: " << 
                     system.variable_name( dbc.vid ) << "("<< dbc.vid << ")" << "=" << dbc.val << 
                     " @ " << bi.get_sideset_name( dbc.bid ) << "(" << dbc.bid << ")";

    ConstFunction<> cf(dbc.val);
    DirichletBoundary bound({dbc.bid}, {dbc.vid}, cf, LOCAL_VARIABLE_ORDER);
    dof_map.add_dirichlet_boundary( bound );
  }

  system.reinit_constraints();
}

/**
 *   Constrain unassinged scalars to avoid singularities.
 */
void ViscoplasticSolver::set_unassigned_scalars() 
{
  DofMap & dof_map = system.get_dof_map();

  // Capture all the scalar vids that have an assigned rigid constrain
  set<uint> assigned;
  for ( auto & sc : curr_bc.scalar ) assigned.insert( sc.svid );

  // Capture all scalars that have not been assigned
  set<uint> unassigned_scalars;
  for ( auto & sv : MODEL->boundary_config.scalars ) {
    const uint vid = system.variable_number (sv.name);
    if ( ! assigned.count( vid ) ) unassigned_scalars.insert(vid);
  }

  vector<dof_id_type> dofi;
  for ( uint vid : unassigned_scalars ) 
  {
    dlog(1) << "Unassigned VID: " << vid;
    dof_map.SCALAR_dof_indices (dofi, vid);
    for ( uint di : dofi ) 
    {
      dlog(1) << "Added constraint row (di=" << di << ").";
      DofConstraintRow cr;
      dof_map.add_constraint_row(di, cr, 0, true);
    }
  }

  
  { // Message 
    ostringstream os; 
    for ( auto vid : unassigned_scalars ) os << system.variable_name(vid) << "(" << vid << ")"; 
    dlog(1) << "Unassigned scalars will be forced to 0: " << os.str();
  }
}

/**
 *
 */
void ViscoplasticSolver::set_scalar_bcs() 
{
  SCOPELOG(1);
  set_unassigned_scalars();

  const MeshBase & mesh = get_mesh();
  const BoundaryInfo & bi = mesh.get_boundary_info();
  DofMap & dof_map = system.get_dof_map();

  // Build a map by bid
  map< sint, vector< BC::ScalarItem > > rigid_map;
  for ( auto & sc : curr_bc.scalar )
  {
    dlog(1) << "Assigning variable '" << system.variable_name(sc.vid) << "' (" << sc.vid << ") to RIGID SCALAR '" << system.variable_name(sc.svid) << "' (" << sc.svid << ").";
    rigid_map[sc.bid].push_back(sc);
  }


  map< uint, uint > n_dofs_by_rvid; /// number of dofs assigned to a scalar (rvid)
  map< uint, set<suint> > bids_by_rvid; /// number of dofs assigned to a scalar (rvid)
  set<uint> dofs_set;
  for (const auto & elem : mesh.active_element_ptr_range()) 
  {
    for (auto side : elem->side_index_range()) 
    {
      vector<sint> vec;

      bi.boundary_ids(elem, side, vec );
      for ( auto bid : vec ) 
      if ( rigid_map.count(bid) ) 
      {
        auto & scalar_vec = rigid_map.at(bid);
        for ( auto & sc : scalar_vec ) {
          auto vid = sc.vid;
          auto rvid = sc.svid;

          vector<dof_id_type> dofi, rdofi;
          dof_map.dof_indices (elem, dofi, vid);
          dof_map.dof_indices (elem, rdofi, rvid);

          bids_by_rvid[rvid].insert(bid); // Reporting info

          for ( uint c : elem->nodes_on_side( side ) ) 
          {
            // Prevents adding a constrain twice in the samce dof
            if ( dofs_set.count(dofi[c]) ) continue;

            // protects lower order variables, which uses only the first nodes of the element.
            if ( c >= dofi.size() ) continue; 

            dofs_set.insert(dofi[c]);

            using util::operator<<;
            DofConstraintRow cr;
            cr.insert( make_pair(rdofi[0], 1) );
            dlog(2) << "Adding constraint row: vid:" << vid << " rvid:" << rvid << " dofi:" << dofi << " rdofi:" << rdofi << "";
            dof_map.add_constraint_row(dofi[c], cr, 0, true);
            
            n_dofs_by_rvid[rvid]++; // Some reporting information
          }
        }
      }
    }
  }

  //
  dlog(1) << "Report on rigid scalar assignment:";
  for ( auto & [ rvid, n ] : n_dofs_by_rvid )
    dlog(1) << "    - " << n << " dofs were assigned to '" << system.variable_name(rvid) << "' (" << rvid << ") ";
  for ( auto & [ rvid, bids ] : bids_by_rvid )
  {
    stringstream os; for ( auto & bid : bids ) os << bi.get_sideset_name(bid) << "(" << bid << ")" << " ";
    dlog(1) << "    - Boundaries assigned to '" << system.variable_name(rvid) << "' (" << rvid << "): " << os.str();
  }
  //
}

/**
 *   Implements the workflow of the system solution.
 *         1. Update dirichlet and rigid constraints
 *         2. call system.solve
 *         3. run sys.get_dof_map().enforce_constraints_exactly(sys) ; sys.update()
 *         4. report convergence
 */
void ViscoplasticSolver::solve()
{
  SCOPELOG(1);
  // Shall we update BCs? 
  if ( curr_bc.update( ts ) ) set_dirichlet_bcs();
  
  es.reinit(); // Maybe only needed if the curr_bc was updated?

  // We might need to rewind in case of TS cuts.
  old_sol = system.solution->clone();

  /** ** ** ** **/
  system.solve();
  /** ** ** ** **/

  bool success = update_adaptive_timestep();
  if ( ! success ) return;  // Do not postprocess if the TS failed

  // Deal with heterogeneous dirichlet boundary constraints
  system.get_dof_map().enforce_constraints_exactly(system);
  system.update();

  ilog1 << "System solved at nonlinear iteration " << system.n_nonlinear_iterations()
    << " , final nonlinear residual norm: " << system.final_nonlinear_residual();

  /** Project the stresses **/
  posproc_stresses();

  /** Export information of the timestep **/
  export_results();

  /** Update for the next timestep **/
  *system.old_local_solution = *system.current_local_solution;
}

/**
 *
 */
void ViscoplasticSolver::do_ts_cut()
{
  SCOPELOG(1);
  ts.cut();

  // Rewind the plastic strain memory to the beginning of the TS
  // TODO: move this TS cutting process to a new method
  MeshBase & mesh = get_mesh();
  for ( const auto & elem : mesh.active_local_element_ptr_range() )
  {
    ViscoPlasticMaterial * mat = get_material( *elem );
    mat->rewind( *elem );
  }

  // Rewind the solution
  *system.solution = *old_sol;
}

/**
 *    If the timestep was cut, return false. Otherwise true;
 */
bool ViscoplasticSolver::update_adaptive_timestep()
{
  if ( system.nonlinear_solver.get() ) 
  { 
    auto * solver = dynamic_cast<PetscNonlinearSolver<Number> *> (system.nonlinear_solver.get());

    if ( solver ) 
    {
      SNESConvergedReason result = solver->get_converged_reason();

      // Converged ?
      if ( result > 0 ) 
      {
        ilog1 << "Nonlinear solver converged (" << SNESConvergedReasons[result] << ").";
        return true;
      }

      // Did not converge. Do the ts cutr.
      wlog1 << "Nonlinear solver diverged (" << SNESConvergedReasons[result] << "). Cutting timestep.";
      do_ts_cut();

      return false;
    } 
  }

  flog << "Could not get the nonlinear solver? Something is terribly wrong.";
  return true;
}

/**
 *
 */
void ViscoplasticSolver::residual_and_jacobian (const NumericVector<Number> & soln, NumericVector<Number> * residual,
                            SparseMatrix<Number> * jacobian, NonlinearImplicitSystem & sys)
{
  Stopwatch sw("ViscoplasticSolver::residual_and_jacobian (ts="+to_string(ts.t_step)+")");
  const DofMap & dof_map = sys.get_dof_map();
  ilog1 << "Assembling poroelastic problem - total DoFs=" << dof_map.n_dofs();

  MeshBase & mesh = get_mesh();

  for ( const auto & elem : mesh.active_local_element_ptr_range() )
  {
    ViscoPlasticMaterial * mat = get_material( *elem );
    mat->residual_and_jacobian( *elem, soln, jacobian, residual );
  }

  harpy_sync_check();

  // Add STOT Boundary Conditions 
  for ( auto & [ elemside, stotitem ] : curr_bc.stot )
  {
    Elem & elem = mesh.elem_ref(elemside.eid);
    // Only on the current processor.
    if ( ! elem.active() ) continue;
    if ( elem.processor_id() != mesh.processor_id() ) continue;

    ViscoPlasticMaterial * mat = get_material( elem );
    ViscoPlasticMaterialBC * bcmat = mat->get_bc_material();
    bcmat->set_bc( stotitem->val );  /// TODO: This is not good.
    bcmat->residual_and_jacobian( elem, elemside.side, soln, jacobian, residual );
  }
  harpy_sync_check();
}

/**
 *
 */
void ViscoplasticSolver::posproc_stresses()
{
  SCOPELOG(1);
  /// Project the stresses into the stress system
  MeshBase & mesh = get_mesh();
  for ( const auto & elem : mesh.active_local_element_ptr_range() )
  {
    ViscoPlasticMaterial * mat = get_material( *elem );
    mat->project_stress( *elem );
  }
  stress_system.solution->close();
}

/**
 *
 */
void ViscoplasticSolver::export_results()
{
  SCOPELOG(1);
  report.do_export();
}

/**
 * Helper stuff 
 */

/** */
ViscoPlasticMaterial * ViscoplasticSolver::get_material( const Elem & elem )
{ uint sid = elem.subdomain_id(); return get_material( sid ); }
/** */
ViscoPlasticMaterial * ViscoplasticSolver::get_material( uint sid )
{
  if  ( ! material_by_sid.count( sid ) ) 
  {
    string sname = get_mesh().subdomain_name( sid );
    flog << "Cannot find material for SID '" << sname << "' (" << sid << ")";
  }
  return material_by_sid.at(sid);
}

}} // ns
