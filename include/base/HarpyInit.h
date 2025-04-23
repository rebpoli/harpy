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

//// verificações em tempo de execução, mesmo em modo opt
#define _harpy_assert_msg(asserted, msg)                               \
  do {                                                                  \
    if (!(asserted)) {                                                  \
      dlog(1) << "[harpy_assert] Assertion `" #asserted "' failed."; \
      libmesh_error_msg(msg);                                           \
    } } while (0)
#define _harpy_assert(asserted) _harpy_assert_msg(asserted, "")

#ifdef DEBUG
#define _harpy_sync_check(comm_obj) do {                            \
    _harpy_assert((comm_obj).verify(std::string(__FILE__).size()));    \
    _harpy_assert((comm_obj).verify(std::string(__FILE__)));           \
    _harpy_assert((comm_obj).verify(__LINE__)); } while (0)
#else
#define _harpy_sync_check(comm_obj) do { } while (0)
#endif

#define harpy_sync_check() _harpy_sync_check((*LIBMESH_COMMUNICATOR))
#define harpy_sync_assert(val) do { \
    _harpy_assert(((*LIBMESH_COMMUNICATOR)).verify(std::string(val)));    \
    } while(0);
