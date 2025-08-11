#pragma once

#include "base/Global.h"

#include "harpy/Solverloop.h"

#include "solver/ViscoplasticSolver.h"

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
    SLViscoplastic( Timestep & ts_ );
    ~SLViscoplastic();

    void solve();

  private:

    Solver * viscoplastic, * thermal, * pressure;
};
