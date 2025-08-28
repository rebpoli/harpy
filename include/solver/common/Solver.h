#pragma once

#include "harpy/Global.h"

#include "libmesh/equation_systems.h"
#include "libmesh/transient_system.h"
#include "libmesh/nonlinear_implicit_system.h"
#include "libmesh/system.h"

#include "libmesh/exodusII_io.h"

#include "util/MpiFileOps.h"

// FWD
namespace libMesh { class Elem ; class MeshBase; } 
namespace timeloop { class Timestep; }
namespace config { class SolverConfig; class BCConfig; }
namespace solver { namespace viscoplastic { class ViscoplasticSolver; } }

/** **/
namespace solver {
namespace common {

/**
 *
 *
 */

class ExplicitMaterial;

using namespace libMesh;

using timeloop::Timestep;
using config::SolverConfig;
using config::BCConfig;
using solver::viscoplastic::ViscoplasticSolver;

/**
 *
 */
class Solver
{
  public:
    Solver( string name_, Timestep & ts_ );
    Solver( Solver * ref, string name_ );
    virtual ~Solver();

    inline MeshBase & get_mesh() { return es.get_mesh(); }
    inline const MeshBase & get_mesh() const { return es.get_mesh(); }

    virtual void init()
      { flog << "Must be defined in the child class."; }
    virtual void solve()
      { flog << "Must be defined in the child class."; }

    void export_exo( string fn );

    string name;
    Timestep & ts;
    bool own_es;

    SolverConfig * config;
    BCConfig & bc_config;

    EquationSystems & es;
};

} } // ns
