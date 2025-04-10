#pragma once

#include "base/Global.h" 

#include "libmesh/fe.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/dof_map.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/dense_submatrix.h"
#include "libmesh/dense_subvector.h"

/**
 *
 * This is an abstract class to be the common inteface 
 * for different materials being passed around in the code.
 *
 * This class holds the FEM shape functions and Gauss quadrature objects.
 *
 * It must be able to build the jacobian and residual within the element.
 *
 */

class MaterialConfig;
class SolverConfig;
class ElemCoupler;
class Coupler;
class Solver;

namespace libMesh { class MeshBase; class System; class Elem; }

using namespace libMesh;

class Material 
{
  public:
    Material( suint sid_, const MaterialConfig & config_ );
    virtual ~Material() {};


    /**
     *   Creates a Material for the subdomain.
     *
     */
    static Material * Factory( suint sid, const MeshBase & mesh, 
                               System & system, const Solver & solver );

    // Common interface to the couplers
    void get_from_element_coupler( string vname, vector<double> & curr , vector<double> & old );
    void get_from_element_coupler( string vname, vector<double> & curr );

    void init_coupler( Elem * elem, ElemCoupler & ec );

    // Interface to any material
    /// Returns a material with the BC definitions and tools
    virtual Material * get_bc_material() 
              { flog << "Must be redifined in the child classes."; return 0; }
    virtual void init_fem() 
              { flog << "Must be redifined in the child classes."; }
    virtual void reinit( const Elem & elem, uint side=255 )
              { flog << "Must be redifined in the child classes."; }
    virtual void reinit( const NumericVector<Number> & soln, const Coupler & coupler, const Elem & elem, uint side=255 )
              { flog << "Must be redifined in the child classes."; }
    virtual void jacobian (const NumericVector<Number> & soln, SparseMatrix<Number> & jacobian )
              { flog << "Must be redifined in the child classes."; }
    virtual void residual (const NumericVector<Number> & soln, NumericVector<Number> & residual )
              { flog << "Must be redifined in the child classes."; }
    virtual bool is_bc() 
              { flog << "Must be redifined in the child classes."; return 0; }
    // Many types of BC can be set. These functions provide a standard interface to pass BCs to any material.
    virtual void set_bc( const RealTensor & value )
              { flog << "Must be redifined in the child classes.";  }

    // Shape functions, quadratures etc
    QGauss qrule;
    vector<dof_id_type> dof_indices;
    unique_ptr<FEBase> fe;  /// The finite element object to hold shape funtions, jxw, etc

    const MaterialConfig & config;
    string name;

  protected:
    void _setup_fem();

    suint sid;   /// Subdomain id

    DenseMatrix<Number> Ke; /// Jacobian for the element
    DenseVector<Number> Re; /// RHS vector for the element

    const ElemCoupler * elem_coupler;

  public:
    /// Must be initialized in the constructor of the child class
    vector< string > required_material_properties; 
};


