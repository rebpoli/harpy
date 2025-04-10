#pragma once

#include "base/Global.h"
#include "harpy/Coupler.h"

#include "libmesh/mesh.h"
#include "libmesh/equation_systems.h"

class Timestep;
class SolverConfig;

/**
 *
 *
 */

class Material;

namespace libMesh { class Elem ; } 
using namespace libMesh;


/**
 *
 */
class Solver
{
  public:
    Solver( string name_, const Timestep & ts_ );

    virtual void solve()
      { flog << "Must be defined in the child class."; }
    virtual void export_results( string basename )
      { flog << "Must be defined in the child class."; }
    virtual MeshBase * get_mesh()
      { flog << "Must be defined in the child class."; return 0; }

    // Adds the products of the current solver to the Coupler,
    // using the target solver material gauss points
    virtual void update_coupler( Coupler & target )
      { flog << "Must be defined in the child class."; }
    virtual void init_trg_coupler( Solver & trg_solver )
      { flog << "Must be defined in the child class."; }
  
    // Initializes the coupler of this object from the material config
    void init_coupler();

    // The material is the same across every solver.
    // Each solver gets its chunk of information as needed
    Material * get_material( const Elem & elem );

    // Directly accessible member. Be careful.
    Coupler coupler;

  protected:
    string name;
    const Timestep & ts;
    map< uint, Material * > material_by_sid;

  public:
    SolverConfig * config;
    Mesh mesh;
    EquationSystems es;

  protected:
};
