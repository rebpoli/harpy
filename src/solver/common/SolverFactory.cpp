
#include "config/ModelConfig.h"
#include "solver/common/SolverFactory.h"
#include "solver/thermal/TSolverCnt.h"
#include "solver/thermal/TSolverFile.h"
#include "solver/pressure/PSolverFile.h"

namespace solver {
namespace common {

using config::MODEL;

/** **/
Solver * SolverFactory::new_thermal( ViscoplasticSolver * ref )
{
  using namespace solver::thermal;

  string sname = "thermal";
  SolverConfig * config = MODEL->solver_config( sname );

  if ( config->external_file.is_defined() ) 
    return new ThermalSolverFromFile( ref,  sname );
  
  return new ThermalSolverConstant( ref,  sname );
}

/** **/
Solver * SolverFactory::new_pressure( ViscoplasticSolver * ref )
{
  using namespace solver::pressure;

  string sname = "pressure";
  SolverConfig * config = MODEL->solver_config( sname );

  if ( config->external_file.is_defined() ) 
    return new PressureSolverFromFile( ref,  sname );
  
  flog << "Unknown pressure solver!";
  return new PressureSolverFromFile( ref,  sname );  // Dummy line, to keep consistency
//  return new PressureSolverConstant( ref,  sname );
}



}} // ns
