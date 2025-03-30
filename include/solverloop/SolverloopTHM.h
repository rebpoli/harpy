#pragma once

#include "base/Global.h"
#include "harpy/Timestep.h"

#include "solver/SolverTHM.h"
#include "harpy/Solverloop.h"

#include "libmesh/equation_systems.h"
#include "libmesh/mesh.h"

/**
 *
 * This is an abstract class.
 *
 */

using namespace libMesh;

class SolverloopTHM : public Solverloop {

  public:
    SolverloopTHM( MeshBase & mesh, const Timestep & ts_ );

    void solve();
    void export_results();

  private:
    EquationSystems es;         // Owned
    const Timestep & ts;        // Owned by Timeloop

    SolverTHM solver;
};
