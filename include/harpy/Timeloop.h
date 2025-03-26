#pragma once

#include "base/Global.h"
#include "base/HarpyInit.h"
#include "harpy/Timestep.h"

#include  "libmesh/mesh.h"

using namespace libMesh;

class Timeloop {
  public:
    Timeloop( MeshBase & mesh_ );
    ~Timeloop();
    void main_loop();

  private:
    MeshBase & mesh;
    Timestep ts;
};

