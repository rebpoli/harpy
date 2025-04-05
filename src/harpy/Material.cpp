
#include "harpy/Material.h"
#include "config/ModelConfig.h" // MODEL global var
#include "config/SolverConfig.h"

#include "libmesh/mesh.h"

/**
 *
 *
 */
Material::Material( const MaterialConfig & mat_conf )
{

}

/**
 *
 */
void Material::reinit()
{
  SCOPELOG(1);
}

/**
 *   Creates a Material for the subdomain.
 *
 */
Material * Material::Factory( uint sid, const MeshBase & mesh, const SolverConfig & svr_config )
{
  SCOPELOG(1);

  set<MaterialConfig> & materials = MODEL->materials;

  string sname = mesh.subdomain_name( sid );

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
  Material * ret = new Material( mat_conf );
  return ret;
}
