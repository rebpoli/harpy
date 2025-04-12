
#include "solver/SolverStress.h"
#include "config/ModelConfig.h"
#include "config/MaterialConfig.h"
#include "material/StressMaterial.h"
#include "harpy/Coupler.h"
#include "libmesh/elem.h"

/**
 *
 */
SolverStress::SolverStress( Solver & ref , string name_ ) :
    Solver( ref, name_ ), system(es.add_system<ExplicitSystem>(name_))
{ init_materials(); }

/**
 *   Creates all the needed materials for the solution (one per subdomain ID).
 *   Initializes the material_by_sid structure.
 */
void SolverStress::init_materials()
{
  SCOPELOG(1);
  MeshBase & mesh = get_mesh();

  set<MaterialConfig> & materials = MODEL->materials;

  // ensures creation of all materials to the current mesh (local elems only)
  for ( const auto & elem : mesh.active_local_element_ptr_range() )
  {
    suint sid = elem->subdomain_id();
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
    material_by_sid[sid] = new StressMaterial( sid, mat_conf, system );
  }
}

/**
 *   Builds a simple structuer from the boundary configuration relating the
 *   material to its current temperature.
 *
 *   This will be used in the update_coupler later.
 */
void SolverStress::solve()
{
  MeshBase & mesh = get_mesh();
  for (const auto & elem : mesh.active_local_element_ptr_range()) 
  {
    uint eid = elem->id();
    if ( ! coupler.count(eid) ) flog << "Element '" << eid << "' missing in coupler";
    ElemCoupler & ec = coupler.at( eid );

    const vector<double> & temp_qp = ec.dbl_params["T"];

    Material * mat = get_material( *elem );
    mat->reinit( *elem );

    // Feeds my own coupler (needed?? normally we would just feed internal vars)
    mat->feed_coupler( ec );

    // TODO --- DO THE PROJECTIONS HERE
  }

}

/**
 *   Initialize the coupler for the target solver with the initial
 *   temperature conditions.
 */
void SolverStress::init_trg_coupler( Solver & trg_solver )
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
void SolverStress::update_coupler( Solver & trg_solver )
{
  SCOPELOG(1);

  Coupler & coupler = trg_solver.coupler;

  // Iterate over the elements and materials of the target solver
  MeshBase & mesh = trg_solver.get_mesh();
  for (const auto & elem : mesh.active_element_ptr_range()) 
  {
    if ( elem->processor_id() != mesh.processor_id() ) continue;  /// Only in the active processor. NOTE: if this requires an element search, that would be a collective task! The loops have to be synchronous.

    uint eid = elem->id();
                                                                   ///
    Material * mat = trg_solver.get_material( *elem );
    mat->reinit( *elem );
    string mname = mat->name;

    if ( ! coupler.count(eid) ) coupler.emplace(eid, ElemCoupler(eid));
    ElemCoupler & ec = coupler.at( eid );

    mat->feed_coupler( ec );
  }
}
