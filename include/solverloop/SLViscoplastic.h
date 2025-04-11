#pragma once

#include "base/Global.h"

#include "harpy/Solverloop.h"

#include "solver/SolverViscoplasticTrial.h"
#include "solver/SolverThermalExplicit.h"

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

  private:

    SolverViscoplasticTrial viscoplastic;
    SolverThermalExplicit thermal;
};
