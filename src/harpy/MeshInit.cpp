
#include "harpy/MeshInit.h"

#include "util/Stopwatch.h"
#include "config/ModelConfig.h"


/**  DEPRECATED --- REMOVE THIS FILE!!! **/

/**
 *
 *
 */
MeshInit::MeshInit( MeshBase & mesh_ ) : mesh(mesh_)
{ 
  load_mesh();
}

/**
 * Loads a mesh from file.
 */
void MeshInit::load_mesh() 
{
  /** Rea mesh from file **/
  {
    string mesh_fn = ""; // MODEL->mesh_filename;  
    Stopwatch sw("mesh.read()");
    ilog1 << "Reading mesh '" << mesh_fn << "' ...";
    mesh.read( mesh_fn );
  }

  /** Convert all elements to the second order 
   * (this might take a while in large meshes) **/
  {
    Stopwatch sw("mesh.all_second_order()");
    mesh.all_second_order();
  }
  
  /** Save the initial mesh in a file for debugging **/
  mesh.write("run/exo/premesh.e");

// Debugging
//  dump(mesh);
}

///**
// *
// *
// *
// */
//void MeshConfig::add_frac_elems( MeshBase & mesh ) {
////  mesh.print_info();
//  BoundaryInfo & bi = mesh.get_boundary_info();
//  auto bid_p = bi.get_id_by_name( "frac_p" );
//  set<boundary_id_type> frac_bids;
//  frac_bids.insert(bid_p);

//  uint n_added=0;
//  class ElAdd { public: uint eid; uint side; vector<Point> pts; };
//  vector<ElAdd> newElems;

//  for (const auto elem : mesh.active_element_ptr_range()) {
//    for (auto side : elem->side_index_range()) {
//      vector<boundary_id_type> vec;
//      bi.boundary_ids(elem, side, vec );
//      for ( auto v : vec ) {
//        // Adiciona elemento 2D igual ao side desse
//        if ( frac_bids.count(v) ) {
//          

//          ElAdd elAdd;
//          elAdd.eid = elem->id();
//          elAdd.side = side;

//          newElems.push_back( elAdd );
//        }
//      }
//    }
//  }

//  // Add new elements
//  map<uint, uint> nodemap; // MAPA: NODE_ID => NOVO NODE_ID
//  for ( ElAdd & elAdd : newElems ) {
//    ++n_added;
//    const Elem * elem = mesh.elem_ptr( elAdd.eid );
//    auto side_elem = elem->build_side_ptr( elAdd.side );
//    Elem * added = mesh.add_elem( Elem::build( side_elem->type() ));
//    added->subdomain_id() = SUBDOMAIN::FRACTURE_2D;
//    for ( uint i=0 ; i<side_elem->n_nodes() ; ++i ) {
//      uint nid = side_elem->node_id( i );
//      if ( ! nodemap.count( nid ) ) {
//        Node * nnode = mesh.add_point( side_elem->point(i) );
//        vector<boundary_id_type> vec;
//        bi.boundary_ids( mesh.node_ptr(nid), vec );
//        for ( auto b : vec ) {
//          dlog(1) << "[AFE] Migrando boundary info ... (bid="<< b <<") nid=" << nid << " @ " << *nnode;
//          bi.add_node(nnode, b);
//        }
//        nodemap[nid] = nnode->id();
//      }

//      added->set_node(i) = mesh.node_ptr( nodemap[nid] );
//    }
//  }

//  mesh.prepare_for_use ();
//  dlog(1) << "Adicionados " << n_added << " elementos de fratura";
//}

///**
// *
// *
// *
// */
//void MeshConfig::set_subdomains( MeshBase & mesh ) {

//  set<string> flow, thermal, frame;
//  CFG.strset( "mesh", "flow", flow );
//  CFG.strset( "mesh", "thermal", thermal );
//  CFG.strset( "mesh", "test_frame", frame );

//  // Feed the SD sets that identifies where each variable is going to reside
//  set<subdomain_id_type> sds; 
//  mesh.subdomain_ids(sds);
//  for ( auto sid : sds ) 
//  {
//    string sname = mesh.subdomain_name( sid );
//    if ( flow.count( sname ) )    SD_PRES.insert( sid );
//    if ( thermal.count( sname ) ) SD_TEMP.insert( sid );
//    if ( frame.count( sname ) ) SD_FRAME.insert( sid );

//    // Currently all domains are part of the CONTINUA
//    SD_CONT.insert(sid);
//  }

//  // Fractures always have pressure
//  SD_PRES.insert( SUBDOMAIN::FRACTURE_2D );

//  // Only for the domains of elements added in this class
//  // The domains set by the users should not be modified.
//  mesh.subdomain_name(SUBDOMAIN::FRACTURE_2D) = "FRAC_2D_ELEMS";
//  mesh.subdomain_name(SUBDOMAIN::WELL) = "WELL";
//  mesh.subdomain_name(SUBDOMAIN::COMPLETION) = "COMPLETION";
//  mesh.subdomain_name(SUBDOMAIN::CASING) = "CASING";

//  // Print debugs
//  dlog(1) << "-- MeshConfig -- domains";
//  for ( auto sid : sds ) dlog(1) << mesh.subdomain_name( sid ) << " => " << sid;
//  dlog(1) << "SD_PRES:" << SD_PRES;
//  dlog(1) << "SD_FRAC:" << SD_FRAC;

//}

///**
// *
// *
// *
// */
//void MeshConfig::dump( const MeshBase & mesh ) {
//  SCOPELOG(9);

//  map<subdomain_id_type, uint> conta_elems;
//  for (const auto & elem : mesh.active_element_ptr_range()) {
//    subdomain_id_type sid = elem->subdomain_id();
//    conta_elems[sid]++;
//  }

//  dlog(1) << "DOMINIOS CONHECIDOS NA MALHA:";
//  std::set< subdomain_id_type > sids;
//  mesh.subdomain_ids (sids);
//  for ( auto sid : sids ) {
//    uint ne = 0;
//    if (conta_elems.count(sid) ) ne = conta_elems[sid];
//    dlog(1) << "    ["<< sid <<"] " << mesh.subdomain_name(sid) << " - " << ne << " elementos.";   
//  }
//  
//  uint sid = mesh.get_id_by_name( "Block" );
//  dlog(1) << "IDBYNAME: Block => " << sid;


//  auto & bi = mesh.get_boundary_info();
////  bi.print_info();
//  set<boundary_id_type> bids = bi.get_boundary_ids();
//  for ( auto bid : bids ) {
//    string name = bi.get_sideset_name(bid);
//    string nname = bi.get_nodeset_name(bid);
//    dlog(1) << "BID:  " << bid << "   / ss: " << name << "  /  ns:  " << nname;
//  }

//  // Constrains dos sides
////    for (auto side : elem->side_index_range()) {
////      vector<boundary_id_type> vec;
////      bi.boundary_ids(elem, side, vec );
////      }

//  for (const auto & elem : mesh.active_element_ptr_range()) {
//    // EDGE CONSTRAINT
//    for (uint edge : elem->edge_index_range()) {
//      vector<boundary_id_type> vec;
//      bi.edge_boundary_ids(elem, edge, vec );
//      for ( auto v : vec ) {
//        string bname    = bi.get_sideset_name( v );
//        dlog(1) << "EDGE CONSTRAINT - ELEM/EDGE: " << elem->id() << "/" << edge << " => " << bname << "(" << v << ")";
//      }
//    }

//    /// SIDE CONSTRAINT
//    for (auto side : elem->side_index_range()) {
//      vector<boundary_id_type> vec;
//      bi.boundary_ids(elem, side, vec );
//      for ( auto v : vec ) {
//        string bname    = bi.get_sideset_name( v );
////        dlog(1) << "SIDE CONSTRAINT - ELEM/SIDE: " << elem->id() << "/" << side << " => " << bname << "(" << v << ")";
//      }
//    }

//    for (auto c : elem->node_index_range() ) {
//      const Node * node = elem->node_ptr(c);
//      vector<boundary_id_type> vec;
//      bi.boundary_ids( node, vec );
////      for ( auto v : vec ) 
////       dlog(1) << "NODE CONSTRAINT: " << v << "  NODE " << node->id() << " @ " << *node;
//    }
//  }


//}


///**
// *
// * Converte o nome dos BOUNDARIES nos BIDS
// *
// */
//set<boundary_id_type> MeshConfig::completion_bids( const libMesh::MeshBase & mesh )
//{
//  auto & bi = mesh.get_boundary_info();

//  set<string> completion;
//  CFG.strset( "mesh", "completion", completion );

//  set<boundary_id_type> ret;
//  for ( auto & s : completion ) 
//  {
//    int bid = bi.get_id_by_name(s);
//    if ( bid == libMesh::BoundaryInfo::invalid_id ) flog << "Boundary id '" << s << "' nao existe na malha. Arquivo de configuraçao invalido."; 
//    ret.insert(bid);
//  }
//  return ret;
//}

///**
// *
// * Converte o nome dos BOUNDARIES nos BIDS
// *
// */
//set<boundary_id_type> MeshConfig::well_bids( const libMesh::MeshBase & mesh )
//{
//  auto & bi = mesh.get_boundary_info();

//  set<string> completion, casing;
//  CFG.strset( "mesh", "completion", completion );
//  CFG.strset( "mesh", "casing", completion );
//  set<string> well;
//  well.insert(completion.begin(), completion.end());
//  well.insert(casing.begin(), casing.end());

//  set<boundary_id_type> ret;
//  for ( auto & s : well ) 
//  {
//    int bid = bi.get_id_by_name(s);
//    if ( bid == libMesh::BoundaryInfo::invalid_id ) flog << "Boundary id '" << s << "' nao existe na malha. Arquivo de configuraçao invalido."; 
//    ret.insert(bid);
//  }
//  return ret;
//}


///**
// *
// * Converte o nome dos dominios nos IDs dos dominios (TERMICO)
// *
// */
//set<uint> MeshConfig::thermal_domains( const libMesh::MeshBase & mesh )
//{
//  set<string> thermal;
//  CFG.strset( "mesh", "thermal", thermal );
//  dlog(1) << "LIDO SET - POROELASTIC::THERMAL:";
//  for ( auto s : thermal ) dlog(1) <<"   " << s;
//  mesh.print_info();
//  set<uint> ret;
//  for ( auto & s : thermal ) 
//  {
//    uint sid = mesh.get_id_by_name( s );
//    if ( sid == libMesh::Elem::invalid_subdomain_id ) flog << "Dominio '" << s << "' nao existe na malha. Arquivo de configuraçao invalido."; 
//    ret.insert(sid);
//  }

//  return ret;
//}
///**
// *
// * Converte o nome dos dominios nos IDs dos dominios (FLUXO)
// *
// */
//set<uint> MeshConfig::flow_domains( const libMesh::MeshBase & mesh )
//{
//  set<string> flow;
//  CFG.strset( "mesh", "flow", flow );
//  dlog(1) << "LIDO SET - POROELASTIC::FLOW:";
//  for ( auto s : flow ) dlog(1) <<"   '" << s << "'";
//  mesh.print_info();

//  set<uint> ret;
//  for ( auto & s : flow ) 
//  {
//    uint sid = mesh.get_id_by_name( s );
//    if ( sid == libMesh::Elem::invalid_subdomain_id ) 
//    {
//      set<subdomain_id_type> sds; mesh.subdomain_ids(sds);
//      dlog(1) << "Subdominios conhecidos:";
//      for ( auto sd : sds ) dlog(1) << "    " << sd << ": " << mesh.subdomain_name(sd);
//      flog << "Dominio '" << s << "' nao existe na malha. Arquivo de configuraçao invalido."; 
//    }
//    ret.insert(sid);
//  }

//  return ret;
//}
