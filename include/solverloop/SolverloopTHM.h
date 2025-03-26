#pragma once

#include "base/Global.h"
#include "harpy/Timestep.h"

#include "libmesh/equation_systems.h"
#include "libmesh/mesh.h"

/**
 *
 * This is an abstract class.
 *
 */

using namespace libMesh;

class SolverloopTHM {

  public:
    SolverloopTHM( MeshBase & mesh, const Timestep & ts_ );
    void solve();

  private:
    EquationSystems es;
    const Timestep & ts;
};
