#pragma once

#include "base/Global.h"

#include "solver/SolverViscoplasticTrial.h"
#include "harpy/Solverloop.h"

class Timestep;

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
