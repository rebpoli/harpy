
#include "harpy/Material.h"
#include "harpy/Solver.h"
#include "config/ModelConfig.h" // MODEL global var
#include "config/SolverConfig.h"
#include "material/ViscoPlasticMaterial.h"

#include "libmesh/elem.h"
#include "libmesh/mesh.h"
#include "libmesh/system.h"

/**
 *
 *
 */
Material::Material( suint sid_, const MaterialConfig & config_ ) :
                          config(config_), name(config.name),
                          sid(sid_), qrule(3), elem_coupler(0)
{ }


/*
 * Define the Material Factory
 */
Material * Material::Factory( suint sid, const MeshBase & mesh, 
                           System & system, const Solver & solver )
{
  SCOPELOG(1);

  set<MaterialConfig> & materials = MODEL->materials;

  string sname = mesh.subdomain_name( sid );

  SolverConfig & svr_config = *( solver.config );
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
  Material * ret = new ViscoPlasticMaterial( sid, mat_conf, system );

  return ret;
}

/**
 *
 *    This function needs the FEBase member to be initialized in the Material to fetch the xyz.
 */
void Material::init_coupler( Elem * elem, ElemCoupler & ec )
{
  // Calc xyz
  const std::vector<Point> & xyz = fe->get_xyz();
  fe->reinit( elem );

  // Feed the coupler
  for ( auto & pname : required_material_properties )
  {
    const MaterialConfig & mconf = config;
    config.get_property( ec.dbl_params[pname], pname, xyz );
  }
}

/**
 *
 */
void Material::get_from_element_coupler( string vname, vector<double> & curr,  vector<double> & old )
{
  if ( ! elem_coupler ) flog << "Element coupler not initialized! Something is wrong.";

  if ( ! elem_coupler->dbl_params.count( vname ) ) flog << "Cannot find variable '" << vname << "' in element coupler. Cannot continue." ;
  if ( ! elem_coupler->dbl_params_old.count( vname ) ) flog << "Cannot find variable '" << vname << "' in element coupler (OLD). Cannot continue." ;

  curr = elem_coupler->dbl_params.at( vname );
  old = elem_coupler->dbl_params_old.at( vname );
}

/**
 *
 */
void Material::get_from_element_coupler( string vname, vector<double> & curr )
{
  if ( ! elem_coupler ) flog << "Element coupler not initialized! Something is wrong.";
  if ( ! elem_coupler->dbl_params.count( vname ) ) {
    dlog(1) << *elem_coupler ;
    flog << "Cannot find variable '" << vname << "' in element coupler. Cannot continue." ;
  }

  curr = elem_coupler->dbl_params.at( vname );
}
