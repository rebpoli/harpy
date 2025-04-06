

#include "solver/SolverViscoplasticTrial.h"
#include "base/HarpyInit.h"
#include "config/ModelConfig.h"
#include "harpy/Timestep.h"
#include "util/MeshUtils.h"


#include "libmesh/elem.h"
#include "libmesh/boundary_info.h"
#include "libmesh/nonlinear_implicit_system.h"
#include "libmesh/transient_system.h"
#include "libmesh/dirichlet_boundaries.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/enum_fe_family.h"
#include "libmesh/enum_order.h"

/**
 *  Creates the system and the materials.
 *
 *  This object owns the mesh and the rquation system.
 */
SolverViscoplasticTrial::SolverViscoplasticTrial( string name_, const Timestep & ts_ ) : 
                   Solver( ts_ ), 
                   name(name_),
                   config (MODEL->solver_config( name ) ) , 
                   mesh( *LIBMESH_COMMUNICATOR ),
                   es( mesh ),
                   system( es.add_system<TransientNonlinearImplicitSystem> ( name ) ),
                   curr_bc( system )
{
  dlog(1) << "SolverViscoplasticTrial: " << *config;

  // Init equationsystems flow
  load_mesh();
  init_materials();
  add_scalar_vars();
  es.init();

  // Init FEM of the materials
  for ( auto & [ sid, mat ] : material_by_sid ) mat->init_fem();
}

/**
 *
 */
void SolverViscoplasticTrial::add_scalar_vars()
{
  for ( auto & sv : MODEL->boundary_config.scalars )
  {
    Order order = Utility::string_to_enum<Order>( sv.order ) ;
    FEFamily fe_family = Utility::string_to_enum<FEFamily>( sv.family );
    system.add_variable( sv.name, order, fe_family );
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
      material_by_sid[sid] = Material::Factory( sid, es.get_mesh(), system, *config );
  }
}

/**
 *   Returns the material for a given element.
 *   Fails if not existing.
 */
Material * SolverViscoplasticTrial::get_material( const Elem & elem, bool reinit )
{
  uint sid = elem.subdomain_id();
  // Consistency check
  if  ( ! material_by_sid.count( sid ) ) 
  {
    string sname = mesh.subdomain_name( sid );
    flog << "Cannot find material for SID '" << sname << "' (" << sid << ")";
  }

  Material * mat = material_by_sid.at(sid);
  if ( reinit ) mat->reinit();

  return mat;
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
//  BCConfig & bcc = MODEL->boundary_config;

//  // BID => VARID => VALUE
//  BCMap<double> dbc;
//  bcs->curr().double_bcs( dbc, _vid_set() );
//  bool drained = bcs->curr().drained;

//  DGConfig dg_config;

//  for ( auto const & me : dbc  ) {
//    auto bid = me.first;
//    auto & vec = me.second;
//    for ( auto & er : vec ) { 
//      auto vid = er.first;
//      auto val = er.second;
//      string bname    = bi.get_sideset_name( bid );
//      string vname    = system.variable_name(vid);

//      ConstFunction<> cf( val );
//      DirichletBoundary bound({bid}, {vid}, cf, LOCAL_VARIABLE_ORDER);
//      dlog1(1) << "["<< fmt_i(ts.t_step) << "] Adding dirichlet boundary condition: " << vname << "("<< vid << ")" << "=" << val << " @ " << bname << "(" << bid << ")";
//      dof_map.add_dirichlet_boundary(bound);
//    }
//  }

//  system.reinit_constraints();

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
  }

}

