#pragma once

#include "base/PetscLogger.h"

#include "libmesh/parallel.h"
#include "libmesh/libmesh.h"

class MyLogInit {
  public:
    MyLogInit();
};

class HarpyInit {
  public:
    HarpyInit( int argc, char ** argv );
    libMesh::Parallel::Communicator & comm() { return lm_init->comm(); }
    PetscLogger pl;
    libMesh::LibMeshInit * lm_init;
    ~HarpyInit();
  private:
    string find_home();
};

