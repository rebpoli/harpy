#include "config/ConfigValidate.h"
#include "config/Config.h"

#include "libmesh/mesh.h"
#include "libmesh/boundary_info.h"

#include <fstream>

namespace config 
{

Config CFG;
using namespace rapidjson;


/**
 *
 * Assistente para obter configurações do json
 *
 */
Config::Config() { }

/**
 *
 *
 */
void Config::init() 
{
  // Arquivo de configuração indicado na variável de ambiente CHIMAS_JSON
  std::string json_fn = getenv("CHIMAS_JSON");

  // Fallback para config.json
  if ( json_fn.empty() ) json_fn = "config.json";

  read_json( _rj, json_fn);

  // Valida o que foi lido
  ConfigValidate v(CFG);
}

/**
 *
 * Verifica se o config.json esta correto.
 *
 */
bool Config::validate() 
{
  SCOPELOG1(1);
  return true;
}

/**
 *
 * Retorna o valor (string) de uma variável de ambiente do sistema.
 *
 */
std::string Config::getenv( const std::string & key ) {
  char * val = std::getenv( key.c_str() );
  string ret;
  if ( val != NULL ) ret = std::string(val);
  return ret;
}


/**
 *
 *
 *
 */
_SolverConfig::_SolverConfig() {
  // PETSc manual - 3.3 KSP: Linear System Solvers - Convergence Tests
  // Default:
  // rtol=1e-5 : the decrease of the residual norm relative to the norm of the right hand side
  // atol=1e-50 : the absolute size of the residual norm
  // dtol=1e5 : the relative increase in the residual
  // maxits=1e4 : the maximum number of allowable iterations
  //
  CFG.dbl("solver","rtol",   rtol, PETSC_DEFAULT);
  CFG.dbl("solver","atol",   atol, PETSC_DEFAULT);
  CFG.dbl("solver","dtol",   dtol, PETSC_DEFAULT);
  CFG.dbl("solver","maxits", maxits, PETSC_DEFAULT);

  CFG.str("solver","ksptype", ksptype, "gmres");
  CFG.str("solver","pctype",  pctype,  "lu");
}

} // ns
