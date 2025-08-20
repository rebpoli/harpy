
#include "solver/ThermalSolverConstant.h"
#include "solver/ViscoplasticSolver.h"
#include "config/ModelConfig.h"
#include "config/BCConfig.h"
#include "harpy/Timestep.h"
#include "util/OutputOperators.h"

#include "material/ViscoPlasticMaterial.h"
#include "libmesh/mesh.h"
#include "libmesh/explicit_system.h"

using harpy_string::iequals;

/*
 *
 */
ThermalSolverConstant::ThermalSolverConstant( Solver * ref_solver_, string name_ ) :
    Solver( ref_solver_, name_ ), system(es.add_system<ExplicitSystem>(name_))
{ 
  SCOPELOG(1);

  setup_variables();
  init_materials(); 
}


/**
 *
 */
void ThermalSolverConstant::setup_variables()
{
  dlog(1) << "Adding variable T ...";
  system.add_variable( "T", FIRST, L2_LAGRANGE ); /// L2_LAGRANGE for discontinuous
}

/**
 *   Creates all the needed materials for the solution (one per subdomain ID).
 *   Initializes the material_by_sid structure.
 */
void ThermalSolverConstant::init_materials()
{
  SCOPELOG(1);
  MeshBase & mesh = get_mesh();

//  set<MaterialConfig> & materials = MODEL->materials;

  // ensures creation of all materials to the current mesh (local elems only)
  for ( const auto & elem : mesh.active_element_ptr_range() )
  {
    suint sid = elem->subdomain_id();
    if  ( material_by_sid.count( sid ) ) continue;

    Material * ref_material = ref_solver->get_material( *elem );
    material_by_sid[sid] = new ExplicitThermalMaterial( ref_material, system );
  }
}

/**
 *   Builds a simple structuer from the boundary configuration relating the
 *   material to its current temperature.
 *
 */
void ThermalSolverConstant::solve()
{
  SCOPELOG(1);
  temperature_by_sid.clear();
  initial_temperature_by_sid.clear();

  MeshBase & mesh = get_mesh();


  ///
  double reftime = bc_config.get_reftime( ts.time ) ;
  if ( ! bc_config.entry_by_time.count( reftime ) ) flog << "Inconsistency in BCConfig. All _reftime_ must have an entry in entry_by_time.";
  const BCConfig::TimeEntry & timeentry = bc_config.entry_by_time.at( reftime );
  ///
  for ( auto & [ sid_name, dbcs ] : timeentry.domain_bcs )
  {
    uint sid = mesh.get_id_by_name(sid_name);
    for ( auto & dbc : dbcs )
    {
      if ( ! iequals ( dbc.vname , "T" ) ) continue;
      temperature_by_sid[ sid ] = dbc.value;
      dlog(1) << "Temperature for '" << sid_name << "(" << sid << ")" << "': " << dbc.value; 
    }
  }

  /// Fetch the initial (reference) time, where the equilibrium temperature is set
  const auto & mapit = bc_config.entry_by_time.cbegin();
  if ( mapit == bc_config.entry_by_time.end() ) flog << "No time entry? Cannot continue.";
  if ( mapit->first >=0 ) flog << "The initial timeentry should be negative, indicating thne initial condition!";
  const BCConfig::TimeEntry & init_timeentry = mapit->second;
  ///
  for ( auto & [ sid_name, dbcs ] : init_timeentry.domain_bcs )
  {
    uint sid = mesh.get_id_by_name(sid_name);
    for ( auto & dbc : dbcs )
    {
      if ( ! iequals ( dbc.vname , "T" ) ) continue;
      initial_temperature_by_sid[ sid ] = dbc.value;
      dlog(1) << "Initial temperature for '" << sid_name << "(" << sid << ")" << "': " << dbc.value; 
    }
  }

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
void ThermalSolverConstant::project_to_system()
{
  SCOPELOG(1);
  MeshBase & mesh = get_mesh();
  for (const auto & elem : mesh.active_element_ptr_range()) 
  {
    ExplicitMaterial * mat = get_explicit_material( *elem );
    mat->reinit( *elem );

    suint sid = elem->subdomain_id();

    string mname = mat->name;
    double temperature = temperature_by_sid[ sid ];

    vector<double> temp_qp;
    uint nqp = mat->qrule.n_points();
    for ( uint qp=0 ; qp<nqp ; qp++ )
      temp_qp.push_back(temperature);

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
void ThermalSolverConstant::update_reference_solver()
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
    suint sid = elem->subdomain_id();

    // init the fe in the explicit material (same qrule as the reference element,
    // different shape function - the one used in project_to_system)
    ExplicitMaterial * mat = get_explicit_material( *elem );
    mat->reinit( *elem ); 

    string mname = mat->name;

    // We need to feed the reference interface
    //
    // TODO: Fix this gambiarra. We need to split the VP interface in thermal, elastic etc and
    //       only propagate what we need.
    Material * vpmat_ = ref_solver->get_material( *elem );
    ViscoPlasticMaterial *vpmat = dynamic_cast< ViscoPlasticMaterial *> ( vpmat_ );

    vpmat->reinit( *elem ); 

    auto & vp_ifc = vpmat->vp_ifc;
    auto & props = vp_ifc.by_elem[eid];

    //
    // TODO: This should work, but it seems that in the first iteration vpmat->qrule.n_points=0?
    //
//     Validation - reference and thermal quadratures must be identical
    if ( vpmat->qrule.n_points() != mat->qrule.n_points() )
      wlog << "Reference and thermal quadrature points are different? Mat.nqp=" << mat->qrule.n_points() << " != VPMat->nqp=" << vpmat->qrule.n_points();

    uint nqp = vpmat->qrule.n_points();
    props.resize(nqp); /// Is this correct?

    /**
     * Update the initial temperature
     */
    {
      double initial_temperature = 0;
      if ( ! initial_temperature_by_sid.count( sid ) ) flog << "No initial temperature defined for material '" << mname << "'.";
      initial_temperature = initial_temperature_by_sid[ sid ];

      // Update the props
      for ( uint qp=0 ; qp< nqp ; qp++ )
        props[qp].initial_temperature = initial_temperature;

      /// Updates the probes
      for ( auto & p : vp_ifc.probes_by_pname_by_elem.probes_by_elem( eid ) )
        p->props.initial_temperature = initial_temperature;
    }

    /**
     *    Update the temperature from the system solution
     */
    {
      /* Quadrature points */
      vector<double> vals_qp;
      mat->eval( vals_qp , "T" );
      // debug
      bool dd=0; // for ( double v : vals_qp ) if ( abs(v-400) > 0.1) dd = 1;
      if ( sid == 1 ) dd = 1;
      if ( dd ) dlog(1) << "[" << mesh.subdomain_name(sid) << "] VALS_QP: " << vals_qp;
      for ( uint qp=0 ; qp< nqp ; qp++ )
        props[qp].temperature = vals_qp[qp];

      /* Update probes */
      for ( auto & p : vp_ifc.probes_by_pname_by_elem.probes_by_elem( eid ) )
        p->props.temperature = mat->eval( p->pt, "T" );
    }
  }
}


/**
 *    DUMPERS
 */
ostream& operator<<(ostream& os, const ThermalSolverConstant & m)
{
  os << "ThermalSolverConstant:" << endl;
  os << "       Name:                 " << m.name << endl;
  os << "       temperature_by_sid:   " << m.temperature_by_sid << endl;
  return os;
}
