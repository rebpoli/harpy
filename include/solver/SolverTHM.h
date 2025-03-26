#pragma once

#include "base/Global.h"

/**
 *
 * Solution for NON - linear THM systems.
 * The underlying materials control which variables are
 * added into the system.
 *
 * The solver name controls the set of configurations
 * that are input.
 *
 * Ex: solver_name = 'poroelastic'
 *                 = 'thermal'
 *
 * This is always nonlinear. The materials
 * must be able to supply a jacobian and a residual for the
 * assigned element.
 *
 */

namespace libMesh { class EquationSystems; }

using namespace libMesh;

class SolverTHM 
{
  public:
    SolverTHM( EquationSystems & es );

  private:
    EquationSystems & es;
};
