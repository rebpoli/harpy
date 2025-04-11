
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
{ 
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
