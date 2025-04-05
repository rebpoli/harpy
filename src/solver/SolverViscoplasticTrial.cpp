
#include "solver/SolverViscoplasticTrial.h"

#include "libmesh/elem.h"
#include "libmesh/boundary_info.h"

#include "config/ModelConfig.h"
#include "base/HarpyInit.h"

/**
 *  Creates the system and the materials.
 *
 *  This object owns the mesh and the rquation system.
 */
SolverViscoplasticTrial::SolverViscoplasticTrial( string name_ ) : Solver(), name(name_),
                                       config (MODEL->solver_config( name ) ) , 
                                       mesh( *LIBMESH_COMMUNICATOR ),
                                       es( mesh )
{
  dlog(1) << "SolverViscoplasticTrial: " << *config;
  load_mesh();
  init_materials();
}

/**
 *    Reads and prepares the mesh for usage.
 */
void SolverViscoplasticTrial::load_mesh()
{
  string fn = config->mesh_filename;
  dlog(1) << "Reading mesh '" << fn << "'...";
  mesh.read( fn );
  dump_mesh();

}

/**
 *   Deletes the owned data structure
 */
SolverViscoplasticTrial::~SolverViscoplasticTrial() 
{
  for ( auto & [ sid, mat ] : material_by_sid ) delete( mat );
  material_by_sid.clear();
}

/**
 *
 */
void SolverViscoplasticTrial::init_materials()
{
  SCOPELOG(1);
  // ensures creation of all materials to the current mesh (local elems only)
  MeshBase & mesh = es.get_mesh();
  for ( const auto & elem : mesh.active_local_element_ptr_range() )
    get_material( *elem );
}

/**
 *   Returns the material for a given element.
 *   Creates a material if not existing.
 */
Material * SolverViscoplasticTrial::get_material( const Elem & elem, bool reinit )
{
  SCOPELOG(1);
  UNUSED(reinit);   /// TODO remove

  CIMap<MaterialConfig> & materials = MODEL->materials;
//  dlog(1) << materials;

  /* ** ** */
  uint sid = elem.subdomain_id();
  if  ( ! material_by_sid.count( sid ) ) 
  {
    MeshBase & mesh = es.get_mesh();

    // 1. Find subdomain name
    string sname = mesh.subdomain_name( sid );

    // 2. Find material name and configuration
    if ( ! config->mat_config_by_name.count( sname ) ) flog << "Cannot find material configuration by name for subdomain '" << sname << "'. The model is inconsistent.";
    auto & mat_conf_id = config->mat_config_by_name.at( sname );

    // 3. Find/build material object
    if ( ! materials.count( mat_conf_id.name ) ) flog << "Cannot find material description for '" << sname << "'. The model is inconsistent.";
    MaterialConfig & mat_conf = materials.at( mat_conf_id.name );

    dlog(1) << "Resoved material:" << mat_conf;

    material_by_sid[sid] = new Material( mat_conf );
//    material_by_sid[sid] = Material::Factory( sid );
  }

  return 0;
}

/**
 *
 */
void SolverViscoplasticTrial::solve()
{

}

/**
 *   Dumps the mesh information.
 */
void SolverViscoplasticTrial::dump_mesh() {
  SCOPELOG(1);
  map<subdomain_id_type, uint> conta_elems;
  for (const auto & elem : mesh.active_element_ptr_range()) {
    subdomain_id_type sid = elem->subdomain_id();
    conta_elems[sid]++;
  }

  dlog(1) << "DOMINIOS CONHECIDOS NA MALHA:";
  std::set< subdomain_id_type > sids;
  mesh.subdomain_ids (sids);
  for ( auto sid : sids ) {
    uint ne = 0;
    if (conta_elems.count(sid) ) ne = conta_elems[sid];
    dlog(1) << "    ["<< sid <<"] " << mesh.subdomain_name(sid) << " - " << ne << " elementos.";   
  }
  
  uint sid = mesh.get_id_by_name( "Block" );
  dlog(1) << "IDBYNAME: Block => " << sid;

  auto & bi = mesh.get_boundary_info();
//  bi.print_info();
  set<boundary_id_type> bids = bi.get_boundary_ids();
  for ( auto bid : bids ) {
    string name = bi.get_sideset_name(bid);
    string nname = bi.get_nodeset_name(bid);
    dlog(1) << "BID:  " << bid << "   / ss: " << name << "  /  ns:  " << nname;
  }

  // Constrains dos sides
//    for (auto side : elem->side_index_range()) {
//      vector<boundary_id_type> vec;
//      bi.boundary_ids(elem, side, vec );
//      }

  for (const auto & elem : mesh.active_element_ptr_range()) {
    // EDGE CONSTRAINT
    for (uint edge : elem->edge_index_range()) {
      vector<boundary_id_type> vec;
      bi.edge_boundary_ids(elem, edge, vec );
      for ( auto v : vec ) {
        string bname    = bi.get_sideset_name( v );
        dlog(1) << "EDGE CONSTRAINT - ELEM/EDGE: " << elem->id() << "/" << edge << " => " << bname << "(" << v << ")";
      }
    }

    /// SIDE CONSTRAINT
    for (auto side : elem->side_index_range()) {
      vector<boundary_id_type> vec;
      bi.boundary_ids(elem, side, vec );
      for ( auto v : vec ) {
        string bname    = bi.get_sideset_name( v );
//        dlog(1) << "SIDE CONSTRAINT - ELEM/SIDE: " << elem->id() << "/" << side << " => " << bname << "(" << v << ")";
      }
    }

    for (auto c : elem->node_index_range() ) {
      const Node * node = elem->node_ptr(c);
      vector<boundary_id_type> vec;
      bi.boundary_ids( node, vec );
//      for ( auto v : vec ) 
//       dlog(1) << "NODE CONSTRAINT: " << v << "  NODE " << node->id() << " @ " << *node;
    }
  }


}
