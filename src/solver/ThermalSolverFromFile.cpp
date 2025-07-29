
#include "solver/ThermalSolverFromFile.h"
#include "solver/ViscoplasticSolver.h"
#include "config/ModelConfig.h"
#include "config/BCConfig.h"
#include "harpy/Timestep.h"
#include "harpy/ExplicitMaterial.h"
#include "util/String.h"
#include "util/OutputOperators.h"

#include "material/ViscoPlasticMaterial.h"
#include "libmesh/mesh.h"
#include "libmesh/explicit_system.h"

#include "util/GridFile.h"

using harpy_string::iequals;

/*
 *
 */
ThermalSolverFromFile::ThermalSolverFromFile( Solver * ref_solver_, string name_ ) :
    Solver( ref_solver_, name_ ), system(es.add_system<ExplicitSystem>(name_)),
    grid(0)
{ 
  SCOPELOG(1);

  SolverConfig * config = MODEL->solver_config(name_);


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
  if ( grid ) delete(grid);
}

/**
 *
 *
 */
void ThermalSolverFromFile::read_from_file()
{
  // Check if we have a complete description from a file
  config->external_file.check(); 

  // Resolve stuff
  string filename = *(config->external_file.filename);
  libMesh::Point origin = *(config->external_file.grid_origin);
  auto grid_type = *(config->external_file.grid_type);

  using enum GRID_TYPE;
  if ( grid_type != RADIAL ) flog << "We only support radial grids so far ...";

  grid = new GridRadialFile(filename);
}

/**
 *   Creates all the needed materials for the solution (one per subdomain ID).
 *   Initializes the material_by_sid structure.
 */
void ThermalSolverFromFile::init_materials()
{
  SCOPELOG(1);
  MeshBase & mesh = get_mesh();

  set<MaterialConfig> & materials = MODEL->materials;

  // ensures creation of all materials to the current mesh (local elems only)
  for ( const auto & elem : mesh.active_element_ptr_range() )
  {
    suint sid = elem->subdomain_id();
    if  ( material_by_sid.count( sid ) ) continue;

    // Get the reference material to copy the quadrature and link
    // to the thermal material
    Material * ref_material = ref_solver->get_material( *elem );

    material_by_sid[sid] = new ExplicitMaterial( ref_material, system );
  }
}

/**
 *
 */
void ThermalSolverFromFile::setup_variables()
{
  dlog(1) << "Adding variable T ...";
  system.add_variable( "T", FIRST, L2_LAGRANGE ); /// L2_LAGRANGE for discontinuous
}

/**
 *   Interpolate the temperature from the external file and assign to the
 *   proper interface.
 *
 */
void ThermalSolverFromFile::solve()
{
  SCOPELOG(1);
//  temperature_by_sid.clear();
//  initial_temperature_by_sid.clear();

//  MeshBase & mesh = get_mesh();


//  ///
//  double reftime = bc_config.get_reftime( ts.time ) ;
//  if ( ! bc_config.entry_by_time.count( reftime ) ) flog << "Inconsistency in BCConfig. All _reftime_ must have an entry in entry_by_time.";
//  const BCConfig::TimeEntry & timeentry = bc_config.entry_by_time.at( reftime );
//  ///
//  for ( auto & [ sid_name, dbcs ] : timeentry.domain_bcs )
//  {
//    uint sid = mesh.get_id_by_name(sid_name);
//    for ( auto & dbc : dbcs )
//    {
//      if ( ! iequals ( dbc.vname , "T" ) ) continue;
//      temperature_by_sid[ sid ] = dbc.value;
//      dlog(1) << "Temperature for '" << sid_name << "(" << sid << ")" << "': " << dbc.value; 
//    }
//  }

//  /// Fetch the initial (reference) time, where the equilibrium temperature is set
//  const auto & mapit = bc_config.entry_by_time.cbegin();
//  if ( mapit == bc_config.entry_by_time.end() ) flog << "No time entry? Cannot continue.";
//  if ( mapit->first >=0 ) flog << "The initial timeentry should be negative, indicating thne initial condition!";
//  const BCConfig::TimeEntry & init_timeentry = mapit->second;
//  ///
//  for ( auto & [ sid_name, dbcs ] : init_timeentry.domain_bcs )
//  {
//    uint sid = mesh.get_id_by_name(sid_name);
//    for ( auto & dbc : dbcs )
//    {
//      if ( ! iequals ( dbc.vname , "T" ) ) continue;
//      initial_temperature_by_sid[ sid ] = dbc.value;
//      dlog(1) << "Initial temperature for '" << sid_name << "(" << sid << ")" << "': " << dbc.value; 
//    }
//  }

//  // Feed the element solver with the temperatures
  project_to_system();
  export_exo("thermal.e");
//  update_reference_solver();
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

  libMesh::Point origin = *(config->external_file.grid_origin);

  MeshBase & mesh = get_mesh();
  for (const auto & elem : mesh.active_element_ptr_range()) 
  {
    ExplicitMaterial * mat = get_explicit_material( *elem );
    mat->reinit( *elem );

    const vector<Point> & xyz_qp = mat->fe->get_xyz();

    string mname = mat->name;


    vector<double> temp_qp;
    uint nqp = mat->qrule.n_points();
    temp_qp.resize(nqp);
    for ( uint qp=0 ; qp<nqp ; qp++ ) 
    {
      Point p = xyz_qp[qp] + origin ;
      temp_qp[qp] = grid->at( ts.time, p(0), p(1), p(2) );
//      dlog(1) << "p:" << p << "  time:" << ts.time << "  T:" << temp_qp[qp];
    } 

    mat->project( temp_qp, "T" );
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

//  /**
//   *   Update the quadrature points
//   */
//  for (const auto & elem : mesh.active_local_element_ptr_range()) 
//  {
//    uint eid = elem->id();
//    suint sid = elem->subdomain_id();

//    // init the fe in the explicit material (same qrule as the reference element,
//    // different shape function - the one used in project_to_system)
//    ExplicitMaterial * mat = get_explicit_material( *elem );
//    mat->reinit( *elem ); 

//    string mname = mat->name;

//    // We need to feed the reference interface
//    //
//    // TODO: Fix this gambiarra. We need to split the VP interface in thermal, elastic etc and
//    //       only propagate what we need.
//    Material * vpmat_ = ref_solver->get_material( *elem );
//    ViscoPlasticMaterial *vpmat = dynamic_cast< ViscoPlasticMaterial *> ( vpmat_ );

//    vpmat->reinit( *elem ); 

//    auto & vp_ifc = vpmat->vp_ifc;
//    auto & props = vp_ifc.by_elem[eid];

//    //
//    // TODO: This should work, but it seems that in the first iteration vpmat->qrule.n_points=0?
//    //
////     Validation - reference and thermal quadratures must be identical
//    if ( vpmat->qrule.n_points() != mat->qrule.n_points() )
//      wlog << "Reference and thermal quadrature points are different? Mat.nqp=" << mat->qrule.n_points() << " != VPMat->nqp=" << vpmat->qrule.n_points();

//    uint nqp = vpmat->qrule.n_points();
//    props.resize(nqp); /// Is this correct?

//    /**
//     * Update the initial temperature
//     */
//    {
//      double initial_temperature = 0;
//      if ( ! initial_temperature_by_sid.count( sid ) ) flog << "No initial temperature defined for material '" << mname << "'.";
//      initial_temperature = initial_temperature_by_sid[ sid ];

//      // Update the props
//      for ( uint qp=0 ; qp< nqp ; qp++ )
//        props[qp].initial_temperature = initial_temperature;

//      /// Updates the probes
//      for ( auto & [ _, m1 ] : vp_ifc.probes_by_pname_by_elem )
//      for ( auto & p : m1[eid] ) 
//        p->props.initial_temperature = initial_temperature;
//    }

//    /**
//     *    Update the temperature from the system solution
//     */
//    {
//      /* Quadrature points */
//      vector<double> vals_qp;
//      mat->eval( vals_qp , "T" );
//      // debug
//      bool dd=0; // for ( double v : vals_qp ) if ( abs(v-400) > 0.1) dd = 1;
//      if ( sid == 1 ) dd = 1;
//      if ( dd ) dlog(1) << "[" << mesh.subdomain_name(sid) << "] VALS_QP: " << vals_qp;
//      for ( uint qp=0 ; qp< nqp ; qp++ )
//        props[qp].temperature = vals_qp[qp];

//      /* Update probes */
//      for ( auto & [ pname, m1 ] : vp_ifc.probes_by_pname_by_elem )
//      for ( auto & p : m1[eid] ) 
//        p->props.temperature = mat->eval( p->pt, "T" );
//    }
//  }
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
