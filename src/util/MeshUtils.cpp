
#include "util/MeshUtils.h"
#include "libmesh/mesh.h"
#include "libmesh/elem.h"
#include "libmesh/boundary_info.h"
using namespace libMesh;

/**
 *   Dumps the mesh information.
 */
void dump_mesh( const MeshBase & mesh ) {
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
