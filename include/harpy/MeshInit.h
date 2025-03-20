#pragma once

#include "base/Global.h"
#include "config/MeshConfig.h"
//#include "libmesh/mesh.h"

namespace libMesh {
  class MeshBase;
}
using namespace libMesh;

class MeshInit
{
  public:
    MeshInit( MeshBase & mesh );
    void load_mesh();
  private:
    MeshBase & mesh;
    MeshConfig config;
};
