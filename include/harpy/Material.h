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
class Solver;

namespace libMesh { class MeshBase; class ExplicitSystem; class Elem; }

using namespace libMesh;

class Material 
{
  public:
    Material( suint sid_, const MaterialConfig & config_ );
    Material( Material & refmat_ );
    virtual ~Material();

    // Interface to any material
    /// Returns a material with the BC definitions and tools
    virtual void init_fem() 
              { flog << "Must be redifined in the child classes."; }

    virtual bool is_bc() { return 0; }  /// Defaults to false. Reimplement in the BC classes to return true;
    virtual string hello() { return "Material."; }
    virtual void update_probes() 
              { flog << "Must be redifined in the child classes."; }

    // Shape functions, quadratures etc
    vector<dof_id_type> dof_indices;
    unique_ptr<FEBase> fe;  /// The finite element object to hold shape funtions, jxw, etc

    Material * refmat;

    const MaterialConfig & config;
    string name;

  protected:
    void _setup_fem();

    suint sid;   /// Subdomain id

    DenseMatrix<Number> Ke; /// Jacobian for the element
    DenseVector<Number> Re; /// RHS vector for the element

  public:

    QGauss qrule;
    const Elem * elem;   /// The element that the material has been reinit'ed to
};

