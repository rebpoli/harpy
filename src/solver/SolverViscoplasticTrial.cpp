

#include "solver/SolverViscoplasticTrial.h"
#include "config/ModelConfig.h"
#include "harpy/Timestep.h"
#include "harpy/DirManager.h"
#include "util/MeshUtils.h"
#include "util/Messages.h"
#include "util/OutputOperators.h"

#include "libmesh/boundary_info.h"
#include "libmesh/elem.h"
#include "libmesh/transient_system.h"
#include "libmesh/nonlinear_implicit_system.h"
#include "libmesh/nonlinear_solver.h"

#include "libmesh/dirichlet_boundaries.h"
#include "libmesh/const_function.h"

#include "libmesh/string_to_enum.h"
#include "libmesh/enum_fe_family.h"
#include "libmesh/enum_order.h"

#include "libmesh/exodusII_io.h"

/**
 *  Creates the system and the materials.
 *
 *  This object owns the mesh and the rquation system.
 */
SolverViscoplasticTrial::SolverViscoplasticTrial( string name_, const Timestep & ts_ ) : 
                   Solver( name_, ts_ ), 
                   system( es.add_system<TransientNonlinearImplicitSystem> ( name ) ),
                   curr_bc( system )
{
  dlog(1) << "SolverViscoplasticTrial: " << *config;

  // Init equationsystems flow
  load_mesh();
  init_materials();
  add_scalar_vars();

  system.nonlinear_solver->residual_object = this;
  system.nonlinear_solver->jacobian_object = this;

  es.init();

  // Init FEM of the materials
  for ( auto & [ sid, mat ] : material_by_sid ) mat->init_fem();
}

/**
 *
 */
void SolverViscoplasticTrial::add_scalar_vars()
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
void SolverViscoplasticTrial::load_mesh()
{
  string fn = config->mesh_filename;
  dlog(1) << "Reading mesh '" << fn << "'...";
  mesh.read( fn );
  dump_mesh( mesh );
}

/**
 *   Deletes the owned data structure
 */
SolverViscoplasticTrial::~SolverViscoplasticTrial() 
{
  for ( auto & [ sid, mat ] : material_by_sid ) delete( mat );
  material_by_sid.clear();
}

/**
 *   Creates all the needed materials for the solution (one per subdomain ID).
 *   Initializes the material_by_sid structure.
 */
void SolverViscoplasticTrial::init_materials()
{
  SCOPELOG(1);
  // ensures creation of all materials to the current mesh (local elems only)
  MeshBase & mesh = es.get_mesh();
  for ( const auto & elem : mesh.active_local_element_ptr_range() )
  {
    uint sid = elem->subdomain_id();
    if  ( ! material_by_sid.count( sid ) ) 
      material_by_sid[sid] = Material::Factory( sid, es.get_mesh(), system, *this );
  }
}

/**
 *
 */
void SolverViscoplasticTrial::set_dirichlet_bcs()
{
  SCOPELOG1(1);

  const MeshBase & mesh = es.get_mesh();
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
void SolverViscoplasticTrial::set_unassigned_scalars() 
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
    for ( uint di : dofi ) {
      DofConstraintRow cr;
      dof_map.add_constraint_row(di, cr, -999, true);
    }
  }

  
  { // Message 
    ostringstream os; 
    for ( auto vid : unassigned_scalars ) os << system.variable_name(vid) << "(" << vid << ")"; 
    dlog(1) << "Unassigned scalars will be forced to -999: " << os.str();
  }
}

/**
 *
 */
void SolverViscoplasticTrial::set_scalar_bcs() 
{
  set_unassigned_scalars();

  const MeshBase & mesh = system.get_mesh();
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

            DofConstraintRow cr;
            cr.insert( make_pair(rdofi[0], 1) );
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
void SolverViscoplasticTrial::solve()
{
  SCOPELOG(1);

  // Shall we update BCs? 
  if ( curr_bc.update( ts ) ) 
  {
    set_dirichlet_bcs();
    set_scalar_bcs();
  }

  es.reinit();
  system.solve();
  system.get_dof_map().enforce_constraints_exactly(system);

  ilog << "System solved at nonlinear iteration " << system.n_nonlinear_iterations()
    << " , final nonlinear residual norm: " << system.final_nonlinear_residual();

  export_exo();
}

/**
 *
 */
void SolverViscoplasticTrial::jacobian 
(const NumericVector<Number> & soln, SparseMatrix<Number> & jacobian, NonlinearImplicitSystem & /*sys*/)
{
  SCOPELOG(1);

  MeshBase & mesh = es.get_mesh();
  for ( const auto & elem : mesh.active_local_element_ptr_range() )
  {
    Material * mat = get_material( *elem );
    mat->reinit( soln, coupler, *elem );
    mat->jacobian( soln, jacobian );
  }

  // Add STOT Boundary Conditions 
  for ( auto & [ elemside, stotitem ] : curr_bc.stot )
  {
    Elem & elem = mesh.elem_ref(elemside.eid);
    // Only on the current processor.
    if ( elem.processor_id() != mesh.processor_id() ) continue;

    Material * mat = get_material( elem );
    Material * bcmat = mat->get_bc_material();
    bcmat->set_bc( stotitem->val );
    bcmat->jacobian( soln, jacobian );
  }
}

/**
 *
 */
void SolverViscoplasticTrial::residual 
(const NumericVector<Number> & soln, NumericVector<Number> & residual, NonlinearImplicitSystem & /*sys*/)
{
  SCOPELOG(1);

  MeshBase & mesh = es.get_mesh();
  for ( const auto & elem : mesh.active_local_element_ptr_range() )
  {
    Material * mat = get_material( *elem );
    mat->reinit( soln, coupler, *elem );
    mat->residual( soln, residual );
  }

  // Add STOT Boundary Conditions 
  for ( auto & [ elemside, stotitem ] : curr_bc.stot )
  {
    Elem & elem = mesh.elem_ref(elemside.eid);
    // Only on the current processor.
    if ( elem.processor_id() != mesh.processor_id() ) continue;

    Material * mat = get_material( elem );
    Material * bcmat = mat->get_bc_material();
    bcmat->reinit( soln, coupler, elem, elemside.side );
    bcmat->set_bc( stotitem->val );  /// TODO: This is not good.
    bcmat->residual( soln, residual );
  }
}

/**
 *   
 *   SUPPORT FUNCTIONS 
 *
 */

/**
 *
 */
void SolverViscoplasticTrial::export_exo()
{
  using namespace harpy_dirmanager;
  string fn = exo_filename( "vptrial", ts );

  ExodusII_IO exo(es.get_mesh());
  exo.write_timestep_discontinuous ( fn, es, 1, ts.time );
}
