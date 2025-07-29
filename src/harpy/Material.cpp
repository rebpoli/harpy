
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
                          refmat(0), config(config_), name(config.name),
                          sid(sid_), qrule(3), elem(0)
{ }

/**
 *  Builds a material inheriting most properties from the reference.
 */
Material::Material( Material * refmat_ ) :
                          refmat( refmat_ ),
                          config(refmat_->config), name(refmat_->name),
                          sid(refmat_->sid), qrule(refmat_->qrule), elem(0)
{ }

/**
 *
 */
Material::~Material() {};


