#pragma once

#include "harpy/Global.h"


namespace libMesh { class MeshBase; }

void dump_mesh( const libMesh::MeshBase & mesh );
