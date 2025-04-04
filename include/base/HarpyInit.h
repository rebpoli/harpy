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

// A global communicator pointer
namespace libMesh { namespace Parallel { class Communicator; } }
extern libMesh::Parallel::Communicator * LIBMESH_COMMUNICATOR;
