
#include "solver/SolverThermalExplicit.h"
#include "config/ModelConfig.h"
#include "config/BCConfig.h"
#include "harpy/Timestep.h"
#include "harpy/Material.h"
#include "material/ThermalMaterial.h"
#include "util/String.h"
#include "util/OutputOperators.h"
#include "solver/ElemProjection.h"

#include "libmesh/mesh.h"
#include "libmesh/explicit_system.h"

using harpy_string::iequals;


/**
 *
 */
SolverThermalExplicit::SolverThermalExplicit( string name_, const Timestep & ts_ ) :
    Solver(name_, ts_ ),
    system(es.add_system<ExplicitSystem>("thermal"))
{

  init_materials();
}

/*
 *
 */
SolverThermalExplicit::SolverThermalExplicit( EquationSystems & es, string name_, const Timestep & ts_ ) :
    Solver(es, name_, ts_ ),
    system(es.add_system<ExplicitSystem>("thermal"))
{
  init_materials();
}

void SolverThermalExplicit::init() 
{
  for ( auto & [ sid, mat ] : material_by_sid ) mat->init_fem();
}

/**
 *   Creates all the needed materials for the solution (one per subdomain ID).
 *   Initializes the material_by_sid structure.
 */
void SolverThermalExplicit::init_materials()
{
  SCOPELOG(1);
  MeshBase & mesh = get_mesh();

  set<MaterialConfig> & materials = MODEL->materials;

  // ensures creation of all materials to the current mesh (local elems only)
  for ( const auto & elem : mesh.active_local_element_ptr_range() )
  {
    suint sid = elem->subdomain_id();
    dlog(1) << "Subdom:" << elem->subdomain_id() << "   /// sid:" << sid;
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

    dlog(1) << "Resoved material:" << mat_conf << " SID:" << sid;
    material_by_sid[sid] = new ThermalMaterial( sid, mat_conf, system );
  }
}

/**
 *   Builds a simple structuer from the boundary configuration relating the
 *   material to its current temperature.
 *
 *   This will be used in the update_coupler later.
 */
void SolverThermalExplicit::solve()
{
  temperature_by_material.clear();
  double reftime = bc_config.get_reftime( ts.time ) ;
  if ( ! bc_config.entry_by_time.count( reftime ) ) flog << "Inconsistency in BCConfig. All _reftime_ must have an entry in entry_by_time.";

  const BCConfig::TimeEntry & timeentry = bc_config.entry_by_time.at( reftime );
  for ( auto & [ material, dbcs ] : timeentry.domain_bcs )
  for ( auto & dbc : dbcs )
  {
    if ( ! iequals ( dbc.vname , "T" ) ) continue;
    temperature_by_material[ material ] = dbc.value;
  }

  // Feed the coupler of the current solver
  update_coupler(*this);

  // Do the projection
  //
  // Explicit solvers do not have materials. It is all done manually.
  //
  MeshBase & mesh = get_mesh();
  for (const auto & elem : mesh.active_local_element_ptr_range()) 
  {
    uint eid = elem->id();
    if ( ! coupler.count(eid) ) flog << "Element '" << eid << "' mission in coupler";
    ElemCoupler & ec = coupler.at( eid );

    const vector<double> & temp_qp = ec.dbl_params["T"];

    Material * mat = get_material( *elem );
    mat->reinit( *elem );
    mat->project( ec, "T" );

  }
  system.solution->close();

  export_exo( "thermal" );
}

/**
 *   Initialize the coupler for the target solver.
 */
void SolverThermalExplicit::init_trg_coupler( Solver & trg_solver )
{
  SCOPELOG(1);

  // For now, the Initial condition is just t=-1
  // So we should be safe by running 2 dry updates so that the we have an identicalk
  // current and old solution for the first timestep
  solve();
  update_coupler( trg_solver ) ;  
  update_coupler( trg_solver ) ;  
}

/**
 *
 */
void SolverThermalExplicit::update_coupler( Solver & trg_solver )
{
  SCOPELOG(1);

  Coupler & coupler = trg_solver.coupler;

  MeshBase & mesh = trg_solver.get_mesh();

  for (const auto & elem : mesh.active_element_ptr_range()) 
  {
    if ( elem->processor_id() != mesh.processor_id() ) continue;  /// Only in the active processor. NOTE: if this requires an element search, that would be a collective task! The loops have to be synchronous.

    uint eid = elem->id();
                                                                   ///
    Material * mat = trg_solver.get_material( *elem );
    mat->reinit( *elem );
    string mname = mat->name;


    // Update temperature
    double temperature = 0;
    if ( ! temperature_by_material.count( mname ) ) wlog << "No temperature defined for material '" << mname << "' at time " << ts.time << " (timestep " << ts.t_step << "). Using 0.";
    temperature = temperature_by_material[ mname ];

    if ( ! coupler.count(eid) ) coupler.emplace(eid, ElemCoupler(eid));

    ElemCoupler & ec = coupler.at( eid );

    // Save old temperature and store the new one
    ec.dbl_params_old["T"] = ec.dbl_params["T"];
    ec.dbl_params["T"].clear();
    for ( uint qp=0 ; qp<mat->qrule.n_points() ; qp++ )
      ec.dbl_params["T"].push_back( temperature );

  }
}

/**
 *    DUMPERS
 */
ostream& operator<<(ostream& os, const SolverThermalExplicit & m)
{
  os << "SolverThermalExplicit:" << endl;
  os << "       Name:                      " << m.name << endl;
  os << "       temperature_by_material:   " << m.temperature_by_material << endl;
  os << "       beta_e_by_material:        " << m.beta_e_by_material << endl;
  os << "       beta_d_by_material:        " << m.beta_d_by_material << endl;
  return os;
}
