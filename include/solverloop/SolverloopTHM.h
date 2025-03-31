#pragma once

#include "base/Global.h"
#include "harpy/Timestep.h"

#include "solver/SolverTHM.h"
#include "harpy/Solverloop.h"

/**
 *
 * This is an abstract class.
 *
 */

using namespace libMesh;

class SolverloopTHM : public Solverloop 
{

  public:
    SolverloopTHM( const Timestep & ts_ );

    void solve();
    void export_results();

  private:

    SolverTHM poroelastic;
};
