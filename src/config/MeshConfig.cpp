#include "config/Config.h"
#include "config/MeshConfig.h"
#include "util/File.h"

#include "libmesh/mesh.h"

using namespace libMesh;

// Instantiate the global stuff
set<libMesh::subdomain_id_type> SD_CONT, SD_PRES, SD_TEMP, SD_FRAME;

/**
 *
 *
 */
MeshConfig::MeshConfig() {
  // Fetch and validate mesh filename
  CFG.str ("mesh","filename",   filename,  "");
  if ( filename.empty() ) flog << "No mesh filename in 'mesh.filename' property.";
  if ( ! file_exists(filename) ) flog << "mesh.filename ('" << filename << "' does not exist.";

  // Fetch any other property as needed.
}

