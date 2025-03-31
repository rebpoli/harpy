
#include "solver/SolverTHM.h"

#include "libmesh/equation_systems.h"

#include "config/ModelConfig.h"

/**
 *  Creates the system and the materials.
 */
SolverTHM::SolverTHM( EquationSystems & es_ ) :
       Solver(), es(es_) 
{
  // Init materials
  map<string, MaterialConfig> & materials = MODEL->materials;
  dlog(1) << materials;
}

/**
 *
 */
void SolverTHM::solve()
{

}
