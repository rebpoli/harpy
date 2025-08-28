
#include "solver/thermal/TSolverFile.h"
#include "solver/viscoplastic/VPSolver.h"
#include "config/ModelConfig.h"
#include "config/BCConfig.h"
#include "timeloop/Timestep.h"
#include "solver/common/ExplicitMaterial.h"

#include "util/String.h"
#include "util/OutputOperators.h"
#include "solver/viscoplastic/VPMaterial.h"

#include "libmesh/mesh.h"
#include "libmesh/elem.h"
#include "libmesh/explicit_system.h"
#include "libmesh/numeric_vector.h"

#include "util/GridFile.h"

namespace solver {
namespace thermal {

using util::iequals;

/*
 *
 */
ThermalSolverFromFile::ThermalSolverFromFile( ViscoplasticSolver * ref_solver_, string name_ ) :
    Solver(  ref_solver_, name_ ),
    ref_solver(ref_solver_),
    system(es.add_system<ExplicitSystem>(name_)),
    grid(0)
{ 
  SCOPELOG(1);

//  SolverConfig * config = MODEL->solver_config(name_);


  setup_variables();

  // Create the postproc materials
  init_materials(); 

  // Read file
  read_from_file();
}

/**
 *
 */
ThermalSolverFromFile::~ThermalSolverFromFile() 
{
  SCOPELOG(1);
  if ( grid ) delete(grid);
  for ( auto & [ sid, mat ] : material_by_sid ) delete( mat );
  material_by_sid.clear();
}

/* */
void ThermalSolverFromFile::init()
{
  for ( auto & [ sid, mat ] : material_by_sid )
    mat->init_fem();
}

/**
 *
 *
 */
void ThermalSolverFromFile::read_from_file()
{
  auto & ef = config->external_file;

  // Check if we have a complete description from a file
  ef.check(); 

  // Resolve stuff
  grid_origin = *(ef.grid_origin);

  string filename = *(ef.filename);
  auto grid_type = *(ef.grid_type);

  using enum GRID_TYPE;
  if ( grid_type != RADIAL ) flog << "We only support radial grids so far ...";

  grid = new GridRadialFile(filename);
  if (ef.min_radius) grid->min_radius = *(ef.min_radius);
}

/**
 *   Creates all the needed materials for the solution (one per subdomain ID).
 *   Initializes the material_by_sid structure.
 */
void ThermalSolverFromFile::init_materials()
{
  SCOPELOG(1);
  MeshBase & mesh = get_mesh();

//  set<MaterialConfig> & materials = MODEL->materials;

  // ensures creation of all materials to the current mesh (local elems only)
  for ( const auto & elem : mesh.active_element_ptr_range() )
  {
    suint sid = elem->subdomain_id();
    if  ( material_by_sid.count( sid ) ) continue;

    // Get the reference material to copy the quadrature and link
    // to the thermal material
    ViscoPlasticMaterial * ref_material = ref_solver->get_material( *elem );

    material_by_sid[sid] = new ExplicitThermalMaterial( ref_material, system );
  }
}

/**
 *
 */
void ThermalSolverFromFile::setup_variables()
{
  dlog(1) << "Adding variable T ...";
  system.add_variable( "T", FIRST, L2_LAGRANGE ); 
  system.add_variable( "T0", FIRST, L2_LAGRANGE ); 
  system.add_variable( "Delta_T", FIRST, L2_LAGRANGE ); 
}

/**
 *   Interpolate the temperature from the external file and assign to the
 *   proper interface.
 *
 */
void ThermalSolverFromFile::solve()
{
  SCOPELOG(1);

  // Feed the element solver with the temperatures
  project_to_system();
  update_reference_solver();
}

/**
 *  Project the temperature to a system. This projection can be
 *  discontinuous (so that we can have steep temperature gradients) or
 *  continuous (where the temperature is going to be smoothen in the
 *  interfaces so that all points are single valued)
 */
void ThermalSolverFromFile::project_to_system()
{
  SCOPELOG(1);

  dlog(1) << "System name: " << system.name();
  dlog(1) << system.get_info();

  MeshBase & mesh = get_mesh();
  for (const auto & elem : mesh.active_element_ptr_range()) 
  {
    ExplicitThermalMaterial * mat = get_material( *elem );
    mat->reinit( *elem );

    const vector<Point> & xyz_qp = mat->fe->get_xyz();
    string mname = mat->name;

    uint nqp = mat->qrule.n_points();

    vector<double> temp_qp, temp0_qp, dtemp_qp;
    temp_qp.resize(nqp);
    temp0_qp.resize(nqp);
    dtemp_qp.resize(nqp);

    for ( uint qp=0 ; qp<nqp ; qp++ ) 
    {
      Point pt = xyz_qp[qp] + grid_origin ;

      double T = grid->at( ts.time, pt(0), pt(1), pt(2) );
      double T0 = grid->at( -1, pt(0), pt(1), pt(2) );
      // This is a hack to avoid noisy temperature
      if ( abs( T - T0 ) < 1 ) T = T0; 

      temp_qp[qp]  = T;
      temp0_qp[qp] = T0;
      dtemp_qp[qp] = T - T0;
    } 

    mat->project( temp_qp,  "T" );
    mat->project( temp0_qp, "T0" );
    mat->project( dtemp_qp, "Delta_T" );
  }

  // Close the system
  system.solution->close(); 
  system.update();
}

/**
 *   Update the reference solver based on the current system solution.
 *
 *   TODO: deal with the initial temperature. The best is likely to have a 
 *         different ThermalSolver for the initial temperature.
 *         
 */
void ThermalSolverFromFile::update_reference_solver()
{
  SCOPELOG(1);

  // The mesh is the same as the reference coupler
  MeshBase & mesh = get_mesh();

  /**
   *   Update the quadrature points
   */
  for (const auto & elem : mesh.active_local_element_ptr_range()) 
  {
    uint eid = elem->id();
//    suint sid = elem->subdomain_id();

    // We need to feed the reference interface
    //
    // TODO: Fix this gambiarra. We need to split the VP interface in thermal, elastic etc and
    //       only propagate what we need.
    ViscoPlasticMaterial * vpmat = ref_solver->get_material( *elem );

    vpmat->reinit( *elem ); 
    const vector<Point> & xyz_qp = vpmat->fe->get_xyz();

    auto & vp_ifc = vpmat->vp_ifc;
    auto & props = vp_ifc.by_elem[eid];

    uint nqp = vpmat->qrule.n_points();
    props.resize(nqp); /// Is this correct?

    /**  Update props **/
    for ( uint qp=0 ; qp< nqp ; qp++ )
    {
      Point pt = xyz_qp[qp] + grid_origin ;
      double T = grid->at( ts.time, pt(0), pt(1), pt(2) );
      double T0 = grid->at( -1, pt(0), pt(1), pt(2) );

      // This is a hack to avoid noisy temperature
      if ( abs( T - T0 ) < 1 ) T = T0; 

      props[qp].temperature = T;
      props[qp].initial_temperature = T0;
    }

    /* Update probes */
    for ( auto & probe : vp_ifc.probes_by_pname_by_elem.probes_by_elem( eid ) )
    {
      Point pt = probe->pt + grid_origin;
      double temp = grid->at( ts.time, pt(0), pt(1), pt(2) );
      double temp0 = grid->at( -1, pt(0), pt(1), pt(2) );
      probe->props.temperature = temp;
      probe->props.initial_temperature = temp0;

//      dlog(1) << "THERMALSOLVERFROMFILE : UPDATE PROBES @ " << pt << "  t=" << ts.time;
//      dlog(1) << "        initial_temperature:" << temp0;
//      dlog(1) << "                temperature:" << temp ;
    }

  }
}


/**
 *    DUMPERS
 */
ostream& operator<<(ostream& os, const ThermalSolverFromFile & m)
{
  os << "ThermalSolverFromFile:" << endl;
  os << "       Name:                 " << m.name << endl;
  return os;
}

/** */
ExplicitThermalMaterial * ThermalSolverFromFile::get_material( const Elem & elem )
{ uint sid = elem.subdomain_id(); return get_material( sid ); }
/** */
ExplicitThermalMaterial * ThermalSolverFromFile::get_material( uint sid )
{
  if  ( ! material_by_sid.count( sid ) ) 
  {
    string sname = get_mesh().subdomain_name( sid );
    flog << "Cannot find material for SID '" << sname << "' (" << sid << ")";
  }
  return material_by_sid.at(sid);
}

}} // ns
