#pragma once

#include "harpy/Global.h"

#include "solverloop/Solverloop.h"

#include "solver/viscoplastic/VPSolver.h"

namespace timeloop { class Timestep; }

/**
 *
 * This is an abstract class.
 *
 */

namespace solverloop {

using solver::viscoplastic::ViscoplasticSolver;
using solver::common::Solver;

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

} // ns
