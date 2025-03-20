#pragma once

#include "base/Global.h"
#include "config/Config.h"

#include "libmesh/id_types.h"

#include <vector>
#include <string>
#include <set>

namespace libMesh {
  class MeshBase;
}

//enum SUBDOMAIN : uint { FRACTURE=900, FRACTURE_2D, DRYROCK, WELL, COMPLETION, CASING, RES, TEST_FRAME, SIDEBURDEN, RES_ISOTHERMAL, DRYROCK_ISOTHERMAL };

// These are only the subdomains that we "create" dynamically.
// All others are set by the user in the mesh
enum SUBDOMAIN : uint { FRACTURE=900, FRACTURE_2D, WELL, COMPLETION, CASING };


// These sets identifies the domains of each variable and are
// set by the user.
extern set<libMesh::subdomain_id_type> SD_CONT, SD_PRES, SD_TEMP, SD_FRAME;

// LISTS OF PROPERTIES OF EACH SUBDOMAIN

// Continuum 
//const set<libMesh::subdomain_id_type> SD_CONT = 
//{
//      SUBDOMAIN::DRYROCK,
//      SUBDOMAIN::DRYROCK_ISOTHERMAL,
//      SUBDOMAIN::RES,
//      SUBDOMAIN::RES_ISOTHERMAL,
//      SUBDOMAIN::SIDEBURDEN,
//      SUBDOMAIN::TEST_FRAME
//};

//// Has pressure (flow)
//const set<libMesh::subdomain_id_type> SD_PRES = 
//{
//      SUBDOMAIN::COMPLETION,
//      SUBDOMAIN::FRACTURE_2D,
//      SUBDOMAIN::RES,
//      SUBDOMAIN::RES_ISOTHERMAL,
//      SUBDOMAIN::SIDEBURDEN,
//      SUBDOMAIN::TEST_FRAME,
//      SUBDOMAIN::WELL
//};

//// Has temperature (non-isothermal)
//const set<libMesh::subdomain_id_type> SD_TEMP = 
//{
//      SUBDOMAIN::CASING,
//      SUBDOMAIN::COMPLETION,
//      SUBDOMAIN::DRYROCK,
//      SUBDOMAIN::FRACTURE_2D,
//      SUBDOMAIN::RES,
//      SUBDOMAIN::TEST_FRAME,
//      SUBDOMAIN::WELL 
//};

// Is fracture
const set<libMesh::subdomain_id_type> SD_FRAC = { SUBDOMAIN::FRACTURE_2D };
// Is well
const set<libMesh::subdomain_id_type> SD_WELL = { SUBDOMAIN::WELL };
// Is comletion
const set<libMesh::subdomain_id_type> SD_COMPLETION = { SUBDOMAIN::COMPLETION };

/**
 *
 *
 */
class MeshConfig 
{
  public:
    MeshConfig();

    string filename;   // The mesh filename to read
};

extern MeshConfig MeshCFG;
