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
    Solver( string name_, const Timestep & ts_ );
    Solver( Solver & ref, string name_ );
    virtual ~Solver();

    inline MeshBase & get_mesh() { return es.get_mesh(); }
    inline const MeshBase & get_mesh() const { return es.get_mesh(); }

    bool same_mesh( const Solver & b ) const { return &get_mesh() == &(b.get_mesh()); }
    void assert_same_mesh( const Solver & b ) const { if ( ! same_mesh(b) ) flog << "Meshes must be equal to continue. TODO: support multimesh."; }

    virtual void solve()
      { flog << "Must be defined in the child class."; }
    virtual void export_results( string basename )
      { UNUSED(basename); flog << "Must be defined in the child class."; }

    virtual void init();

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
};
