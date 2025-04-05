#pragma once

#include "base/Global.h"
#include "harpy/Timestep.h"

#include "solver/SolverViscoplasticTrial.h"
#include "harpy/Solverloop.h"

/**
 *
 * This is an abstract class.
 *
 */

using namespace libMesh;

class SLViscoplastic : public Solverloop 
{

  public:
    SLViscoplastic( const Timestep & ts_ );

    void solve();
    void export_results();

  private:

    SolverViscoplasticTrial viscoplastic;
};
