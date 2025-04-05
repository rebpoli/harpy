#pragma once

#include "base/Global.h"


namespace libMesh { class MeshBase; }

void dump_mesh( const libMesh::MeshBase & mesh );
