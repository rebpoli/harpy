
#include "harpy/Material.h"
#include "harpy/Solver.h"
#include "config/ModelConfig.h" // MODEL global var
#include "config/SolverConfig.h"
#include "material/ViscoPlasticMaterial.h"

#include "libmesh/elem.h"
#include "libmesh/mesh.h"
#include "libmesh/explicit_system.h"

#include <regex>


// Some local Regex
const string prop_name = R"(([-a-zA-Z_0-9]+))";
const regex RE_PROPERTY ( R"(^)" + prop_name + R"(\.)" + prop_name + R"($)");

/**
 *
 *
 */
Material::Material( suint sid_, const MaterialConfig & config_ ) :
                          config(config_), name(config.name),
                          sid(sid_), qrule(3), elem(0)
{ }

/**
 *  Builds a material inheriting most properties from the reference.
 */
Material::Material( Material & refmat ) :
                          config(refmat.config), name(refmat.name+"-derived"),
                          sid(refmat.sid), qrule(refmat.qrule), elem(0)
{ }

/**
 *
 */
Material::~Material() {};


