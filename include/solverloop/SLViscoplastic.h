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

    void load_sig0_file( string filename );
    void solve();

  private:

    ViscoplasticSolver *viscoplastic;
    Solver *thermal, *pressure;

    friend restart::File;
};
