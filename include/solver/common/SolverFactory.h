
#pragma once 

#include "harpy/Global.h"

namespace solver { namespace viscoplastic { class ViscoplasticSolver; } }

namespace solver {
namespace common {

class Solver;
using solver::viscoplastic::ViscoplasticSolver;


/**
 *   Factories of solvers.
 */
struct SolverFactory
{
  static Solver * new_thermal( ViscoplasticSolver * ref );
  static Solver * new_pressure( ViscoplasticSolver * ref );
};


}} // ns
