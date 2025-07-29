#pragma once

#include "base/Global.h"

#include "libmesh/equation_systems.h"
#include "libmesh/transient_system.h"
#include "libmesh/nonlinear_implicit_system.h"
#include "libmesh/system.h"

#include "libmesh/exodusII_io.h"

class Timestep;
class SolverConfig;

/**
 *
 *
 */

class Material;
class ExplicitMaterial;
class BCConfig;

namespace libMesh { class Elem ; class MeshBase; } 
using namespace libMesh;


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

    virtual void solve()
      { flog << "Must be defined in the child class."; }

    /// The material can be provided for the dependent solvers
    Material * get_material( const Elem & elem );
    Material * get_material( uint sid );

    virtual void init();

    void export_exo( string fn );

    string name;
    Timestep & ts;
    map< uint, Material * > material_by_sid;
    bool own_es;

    SolverConfig * config;
    BCConfig & bc_config;

    EquationSystems & es;

    Solver * ref_solver;
};


/**
 *   Factories of solvers.
 */
struct SolverFactory
{
  static Solver * new_thermal( Solver * ref );
};
