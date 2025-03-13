#pragma once

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

namespace geo
{

  // Propriedades de uma fratura
  typedef struct
  {
    // id é um sufixo único para identificar as entidades
    // correspondentes à fratura: id_point, id_line, id_face.
    std::string id;
    // length define o comprimento da fratura
    double length;
    // depth define a profundidade da fratura
    double depth;
    // Altura dos elementos planos acima e abaixo da fratura
    double faceLayer;
    // faceMeshSize define o tamanho aproximado dos elementos
    // na face da fratura
    double faceMeshSize;
    // tipMeshSize define o tamanho aproximado dos elementos
    // na ponta da fratura
    double tipMeshSize;
    // Raio de refino
    double meshFieldRadius;
  } Crack;

}

/**
 *
 *
 */
class MeshConfig 
{
  public:
    MeshConfig();
    vector<geo::Crack> cracks;
    string filename() { return "run/msh/" + CFG.str("mesh", "file"); }
    string model() { return CFG.str("mesh", "model"); }

    void init();
    void load_mesh( libMesh::MeshBase & mesh );
    void add_frac_elems( libMesh::MeshBase & mesh );
    void add_well_elems( libMesh::MeshBase & mesh );
    void set_subdomains( libMesh::MeshBase & mesh );
    void dump( const libMesh::MeshBase & mesh );

    set<libMesh::boundary_id_type> completion_bids( const libMesh::MeshBase & mesh );
    set<libMesh::boundary_id_type> well_bids( const libMesh::MeshBase & mesh );
    set<uint> flow_domains( const libMesh::MeshBase & mesh );
    set<uint> thermal_domains( const libMesh::MeshBase & mesh );

    const map<uint,uint> & get_well_mirrors() const { return well_mirrors; }

  private:
    void _init_cracks();

    /// Associacao entre os elementos (EDGE) de poço e de completação ou casing
    map<uint,uint> well_mirrors;
};

extern MeshConfig MeshCFG;
