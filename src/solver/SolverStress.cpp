
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
  for ( const auto & elem : mesh.active_element_ptr_range() )
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

    dlog(1) << "Resoved Stress material:" << mat_conf << " SID:" << sid;
    material_by_sid[sid] = new StressMaterial( sid, mat_conf, system );
  }
}

/**
 *   Builds a simple structuer from the boundary configuration relating the
 *   material to its current temperature.
 *
 */
void SolverStress::solve()
{
  SCOPELOG(1);
  // Auto-feed our coupler to calculate the stresses
  update_coupler( *this );

  MeshBase & mesh = get_mesh();
  for (const auto & elem : mesh.active_element_ptr_range()) 
  {
    MaterialExplicit * mat = get_explicit_material( *elem );

    mat->reinit( coupler, *elem );
    ElemCoupler & ec = coupler.elem_coupler( elem->id() );
    mat->project_tensor( ec, "sigeff" );
    mat->project_tensor( ec, "sigtot" );
    mat->project_tensor( ec, "deviatoric" );
    mat->project( ec, "von_mises" );
    mat->project_tensor( ec, "plastic_strain" );
  }

}

/**
 *    Feeds the target solver with the stress measurements.
 *
 *    The information flows in this direction THIS_SOLVER => TRG_SOLVER
 *
 */
void SolverStress::update_coupler( Solver & trg_solver )
{
  SCOPELOG(1);
  assert_same_mesh( trg_solver );

  MeshBase & trg_mesh = trg_solver.get_mesh();
  for (const auto & src_elem : trg_mesh.active_element_ptr_range())  // This is a collective loop (no local here please!)
  {
    uint eid = src_elem->id();
    ElemCoupler & trg_ec = trg_solver.coupler.elem_coupler( eid );
    ElemCoupler & src_ec = coupler.elem_coupler( eid );
    Material * src_mat = get_material( *src_elem );
    src_mat->reinit( coupler, *src_elem );
    src_mat->feed_coupler( trg_ec );
  }
}
