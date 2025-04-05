#pragma once

#include "base/Global.h" 

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

namespace libMesh { class MeshBase; }
using namespace libMesh;

class Material 
{
  public:
    Material( const MaterialConfig & mat_conf );

    void reinit();

    static Material * Factory( uint sid, const MeshBase & mesh, const SolverConfig & svr_config );
};
