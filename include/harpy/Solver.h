#pragma once

#include "base/Global.h"
#include "harpy/Coupler.h"

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
class BCConfig;

namespace libMesh { class Elem ; class MeshBase; } 
using namespace libMesh;


/**
 *
 */
class Solver
{
  public:
    Solver( string name_, const Timestep & ts_ );
    Solver( Solver & ref, string name_ );
    virtual ~Solver();

    inline MeshBase & get_mesh() { return es.get_mesh(); }

    virtual void solve()
      { flog << "Must be defined in the child class."; }
    virtual void export_results( string basename )
      { flog << "Must be defined in the child class."; }

    virtual void init();

    // Initializes the coupler of this object from the material config
    void init_coupler();

    // The material is the same across every solver.
    // Each solver gets its chunk of information as needed
    Material * get_material( const Elem & elem );
    virtual void init_materials()
      { flog << "Must be defined in the child class."; }
    void export_exo( string fn );

    string name;
    const Timestep & ts;
    map< uint, Material * > material_by_sid;
    bool own_es;

    SolverConfig * config;
    BCConfig & bc_config;

    EquationSystems & es;

    // Directly accessible member. Be careful.
    Coupler coupler;

};
