#include "base/Global.h"
#include "config/Config.h"
#include "config/MeshConfig.h"
#include "util/Stopwatch.h"
#include "util/OutputOperators.h"

#include "libmesh/mesh.h"
#include "libmesh/elem.h"
#include "libmesh/boundary_info.h"

using namespace libMesh;
using namespace geo;

// Instantiate the global stuff
set<libMesh::subdomain_id_type> SD_CONT, SD_PRES, SD_TEMP, SD_FRAME;

MeshConfig MeshCFG;

/**
 *
 *
 */
MeshConfig::MeshConfig() 
{
}

void MeshConfig::init() 
{
  _init_cracks();
}

/**
 *
 *
 */
void MeshConfig::load_mesh( MeshBase & mesh ) {
  /// Le malha do arquivo
  string fn = MeshConfig().filename(); 
  ilog1 << "Lendo malha '" << fn << "'";
  mesh.read( fn );
  dump(mesh);

  /// Converte elementos para segunda ordem
  {
    Stopwatch sw("mesh.all_second_order()");
//    mesh.all_second_order();
    mesh.all_second_order();
  }
  
  // GAMBIARRA!
  if ( MeshConfig().model() == "poco" ){
    mesh.get_boundary_info().sideset_name(8) = "well0";
    mesh.get_boundary_info().sideset_name(9) = "well1";
    mesh.get_boundary_info().sideset_name(10) = "well2";
  }

  // Gambiarra para conhecer a fenda da fratura
  mesh.get_boundary_info().sideset_name(20000) = "frac_open";

  /// Ajusta subdominios
  set_subdomains(mesh);

  add_well_elems(mesh);
  add_frac_elems(mesh);

  // GAMBIARRA!
  if ( MeshConfig().model() == "poco" ){
    mesh.get_boundary_info().sideset_name(8) = "well0";
    mesh.get_boundary_info().sideset_name(9) = "well1";
    mesh.get_boundary_info().sideset_name(10) = "well2";
  }

  // Gambiarra para conhecer a fenda da fratura
  mesh.get_boundary_info().sideset_name(20000) = "frac_open";


  /// Salva a malha inicial em disco
  mesh.write("run/exo/premesh.e");

  // Depuracao
  dump(mesh);
}

/**
 *
 *
 *
 */
void MeshConfig::add_well_elems( libMesh::MeshBase & mesh )
{
  libMesh::BoundaryInfo & bi = mesh.get_boundary_info();

  // BIDs que definem o poco: os edges marcados com essas tags são poços
  set<boundary_id_type> wellbids = well_bids(mesh);
//  wellbids.insert( bi.get_id_by_name( "well0" ) );
//  wellbids.insert( bi.get_id_by_name( "well1" ) );
//  wellbids.insert( bi.get_id_by_name( "well2" ) );
  dlog(1) << "WELL BIDS:";
  for ( auto bid : wellbids ) dlog(1) << "    " << bid;

  //
  // 1. ADICIONA ELEMENTOS LINHA PARA REPRESETNAR OS INJECTION POINTS
  //
  // 1.1 IDENTIFICA OS ELEMENTOS PARA ADICIONAR
  //
  struct WellE { boundary_id_type bid; ElemType etype; set<int> nids; vector<int> nid_vec; } ;
  struct WellECompare {
    bool operator()( const WellE & a, const WellE & b ) const {
      if ( a.bid < b.bid ) return true;
      if ( a.bid > b.bid ) return false;

      if ( a.nids < b.nids ) return true;
      if ( a.nids > b.nids ) return false;

      return false;
    }
  };
  set< WellE, WellECompare > welle;

  for (const auto & elem : mesh.active_element_ptr_range()) {
    if ( ! SD_CONT.count( elem->subdomain_id() ) ) continue; // analisa os elementos do continuo

//    uint n_edges = elem->n_edges();

    for ( uint eid : elem->edge_index_range() )
    {
      std::unique_ptr<Elem> edge = elem->build_edge_ptr(eid);
      map<boundary_id_type, uint> ncnt; // pelo menos 2 nos em um bid para caracterizar um edge

      for (auto nid : edge->node_index_range()) 
      {
        vector<boundary_id_type> vec;
        bi.boundary_ids( edge->node_ptr(nid), vec );
        for ( auto bid : vec ) ncnt[bid] += 1;
      }

      // verifica se tem mais de um nó com esse bid no edge
      for ( auto p : ncnt ) {
        if ( p.second < 2 ) continue;
        boundary_id_type bid = p.first;

//        dlog(1) << "   EDGE COM 2 NODOS COM BID=" << bid;

        // So adiciona edges de poço
        if ( ! wellbids.count(bid) ) continue;
        
        WellE we = { bid, edge->type(), {}, {} };

        // Aqui sabemos que tem que adicionar um elemento 1D na malha no dominio INJ_POINT!
        dlog(1) << "ENCONTREI '" << p.second << "' NOS COM A BID " << bid << " @ ELEM:" << elem->id() << " / Edge:" << eid;
        for ( uint i=0 ; i<edge->n_nodes() ; ++i ) {
          uint nid = edge->node_id( i );
          we.nids.insert(nid);
          we.nid_vec.push_back(nid);
        }

        // Como aumentamos a ordem do elemento, é importante colocar o BID em todos
        for (auto nid : edge->node_index_range()) 
          bi.add_node( edge->node_ptr(nid), bid );

        // Coloca no conjunto para manter unico
        if ( welle.count(we) ) dlog(1) << "CONTAINER JA TINHA O WE!";
        welle.insert(we);
      }
    }
  }

  /// Aqui tem que ter um procedimento para eliminar edges iguais

  // Novos nós, do elemento de poço gemeo
  map<uint, uint> nodemap; // MAPA: NODE_ID => NOVO NODE_ID
  set<boundary_id_type> comp_bids = completion_bids(mesh);

  //
  // ADICIONA OS NOVOS ELEMENTOS DE POÇO, CASING E COMPLETACAO
  //
  for ( auto we : welle ) 
  {
    // ELEMENTO DE COMPLETAÇÃO OU CASING, COM OS NÓS EXISTENTES
    Elem * compl_elem = mesh.add_elem( Elem::build( we.etype ));
    uint sid = SUBDOMAIN::COMPLETION;
    if ( ! comp_bids.count( we.bid ) ) sid = SUBDOMAIN::CASING;
    compl_elem->subdomain_id() = sid;
    for ( uint i=0 ; i<compl_elem->n_nodes() ; ++i ) {
      uint nid = we.nid_vec[i];
      compl_elem->set_node(i) = mesh.node_ptr( nid );
    }

    // ELEMENTO DE POÇO, COM NÓS NOVOS
    Elem * well_elem = mesh.add_elem( Elem::build( we.etype ));
    well_elem->subdomain_id() = SUBDOMAIN::WELL;
    // Adiciona novos nos
    for ( uint i=0 ; i<well_elem->n_nodes() ; ++i ) {
      uint nid = we.nid_vec[i];
      if ( ! nodemap.count( nid ) ) {
        Node * nnode = mesh.add_point( mesh.point(nid) );
        nodemap[nid] = nnode->id();

        // Migra boundary info ...
        vector<boundary_id_type> vec;
        bi.boundary_ids( mesh.node_ptr(nid), vec );
        for ( auto b : vec ) {
          dlog(1) << "[AWE] Migrando boundary info ... (bid="<< b <<") nid=" << nid << " @ " << *nnode << "  well_elem:" << well_elem->id();
          bi.add_node(nnode, b);
        }

      }
      well_elem->set_node(i) = mesh.node_ptr( nodemap[nid] );
    }

    // Anota a associacao entre os elementos, bidirecional
    well_mirrors[ well_elem->id() ] = compl_elem->id();;
    well_mirrors[ compl_elem->id() ] = well_elem->id();;
  }
}

/**
 *
 *
 *
 */
void MeshConfig::add_frac_elems( MeshBase & mesh ) {
//  mesh.print_info();
  BoundaryInfo & bi = mesh.get_boundary_info();
  auto bid_p = bi.get_id_by_name( "frac_p" );
  set<boundary_id_type> frac_bids;
  frac_bids.insert(bid_p);

  uint n_added=0;
  class ElAdd { public: uint eid; uint side; vector<Point> pts; };
  vector<ElAdd> newElems;

  for (const auto elem : mesh.active_element_ptr_range()) {
    for (auto side : elem->side_index_range()) {
      vector<boundary_id_type> vec;
      bi.boundary_ids(elem, side, vec );
      for ( auto v : vec ) {
        // Adiciona elemento 2D igual ao side desse
        if ( frac_bids.count(v) ) {
          

          ElAdd elAdd;
          elAdd.eid = elem->id();
          elAdd.side = side;

          newElems.push_back( elAdd );
        }
      }
    }
  }

  // Add new elements
  map<uint, uint> nodemap; // MAPA: NODE_ID => NOVO NODE_ID
  for ( ElAdd & elAdd : newElems ) {
    ++n_added;
    const Elem * elem = mesh.elem_ptr( elAdd.eid );
    auto side_elem = elem->build_side_ptr( elAdd.side );
    Elem * added = mesh.add_elem( Elem::build( side_elem->type() ));
    added->subdomain_id() = SUBDOMAIN::FRACTURE_2D;
    for ( uint i=0 ; i<side_elem->n_nodes() ; ++i ) {
      uint nid = side_elem->node_id( i );
      if ( ! nodemap.count( nid ) ) {
        Node * nnode = mesh.add_point( side_elem->point(i) );
        vector<boundary_id_type> vec;
        bi.boundary_ids( mesh.node_ptr(nid), vec );
        for ( auto b : vec ) {
          dlog(1) << "[AFE] Migrando boundary info ... (bid="<< b <<") nid=" << nid << " @ " << *nnode;
          bi.add_node(nnode, b);
        }
        nodemap[nid] = nnode->id();
      }

      added->set_node(i) = mesh.node_ptr( nodemap[nid] );
    }
  }


//          for ( uint i=0 ; i<side_elem->n_nodes() ; ++i ) {
//            // mapeia o node id
//            uint nid = side_elem->node_id( i );

//            // adiciona pontos e ajusta boundary infos
//            if ( ! nodemap.count( nid ) ) {
//              Node * nnode = mesh.add_point( side_elem->point(i) );
//              vector<boundary_id_type> vec;
//              bi.boundary_ids( mesh.node_ptr(nid), vec );
//              for ( auto b : vec ) {
//                dlog(1) << "[AFE] Migrando boundary info ... (bid="<< b <<") nid=" << nid << " @ " << *nnode;
//                bi.add_node(nnode, b);
//              }
//              nodemap[nid] = nnode->id();
//            }

//            added->set_node(i) = mesh.node_ptr( nodemap[nid] );
//          }

//          goto next_elem; // se livrando de 2 loops
//        }
//      }
//    }
//    next_elem: continue;


  mesh.prepare_for_use ();
  dlog(1) << "Adicionados " << n_added << " elementos de fratura";
}

/**
 *
 *
 *
 */
void MeshConfig::set_subdomains( MeshBase & mesh ) {

  set<string> flow, thermal, frame;
  CFG.strset( "mesh", "flow", flow );
  CFG.strset( "mesh", "thermal", thermal );
  CFG.strset( "mesh", "test_frame", frame );

  // Feed the SD sets that identifies where each variable is going to reside
  set<subdomain_id_type> sds; 
  mesh.subdomain_ids(sds);
  for ( auto sid : sds ) 
  {
    string sname = mesh.subdomain_name( sid );
    if ( flow.count( sname ) )    SD_PRES.insert( sid );
    if ( thermal.count( sname ) ) SD_TEMP.insert( sid );
    if ( frame.count( sname ) ) SD_FRAME.insert( sid );

    // Currently all domains are part of the CONTINUA
    SD_CONT.insert(sid);
  }

  // Fractures always have pressure
  SD_PRES.insert( SUBDOMAIN::FRACTURE_2D );

  // Only for the domains of elements added in this class
  // The domains set by the users should not be modified.
  mesh.subdomain_name(SUBDOMAIN::FRACTURE_2D) = "FRAC_2D_ELEMS";
  mesh.subdomain_name(SUBDOMAIN::WELL) = "WELL";
  mesh.subdomain_name(SUBDOMAIN::COMPLETION) = "COMPLETION";
  mesh.subdomain_name(SUBDOMAIN::CASING) = "CASING";

  // Print debugs
  dlog(1) << "-- MeshConfig -- domains";
  for ( auto sid : sds ) dlog(1) << mesh.subdomain_name( sid ) << " => " << sid;
  dlog(1) << "SD_PRES:" << SD_PRES;
  dlog(1) << "SD_FRAC:" << SD_FRAC;

}

/**
 *
 *
 *
 */
void MeshConfig::dump( const MeshBase & mesh ) {
  SCOPELOG(9);

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

/**
 *
 *
 */
void MeshConfig::_init_cracks() {
  // Carregar configurações de fraturas
  if ( CFG.exists("mesh", "cracks" ) ) {
    const rapidjson::Value& crack_config = CFG._rj["mesh"]["cracks"];
    for ( auto& crack : crack_config.GetArray() )
    {

      cracks.push_back(
      {
        crack["id"].GetString(),
        crack["length"].GetDouble(),
        crack["depth"].GetDouble(),
        crack["faceLayer"].GetDouble(),
        crack["faceMeshSize"].GetDouble(),
        crack["tipMeshSize"].GetDouble(),
        crack["meshFieldRadius"].GetDouble()
      });
    }
  }

}


/**
 *
 * Converte o nome dos BOUNDARIES nos BIDS
 *
 */
set<boundary_id_type> MeshConfig::completion_bids( const libMesh::MeshBase & mesh )
{
  auto & bi = mesh.get_boundary_info();

  set<string> completion;
  CFG.strset( "mesh", "completion", completion );

  set<boundary_id_type> ret;
  for ( auto & s : completion ) 
  {
    int bid = bi.get_id_by_name(s);
    if ( bid == libMesh::BoundaryInfo::invalid_id ) flog << "Boundary id '" << s << "' nao existe na malha. Arquivo de configuraçao invalido."; 
    ret.insert(bid);
  }
  return ret;
}

/**
 *
 * Converte o nome dos BOUNDARIES nos BIDS
 *
 */
set<boundary_id_type> MeshConfig::well_bids( const libMesh::MeshBase & mesh )
{
  auto & bi = mesh.get_boundary_info();

  set<string> completion, casing;
  CFG.strset( "mesh", "completion", completion );
  CFG.strset( "mesh", "casing", completion );
  set<string> well;
  well.insert(completion.begin(), completion.end());
  well.insert(casing.begin(), casing.end());

  set<boundary_id_type> ret;
  for ( auto & s : well ) 
  {
    int bid = bi.get_id_by_name(s);
    if ( bid == libMesh::BoundaryInfo::invalid_id ) flog << "Boundary id '" << s << "' nao existe na malha. Arquivo de configuraçao invalido."; 
    ret.insert(bid);
  }
  return ret;
}


/**
 *
 * Converte o nome dos dominios nos IDs dos dominios (TERMICO)
 *
 */
set<uint> MeshConfig::thermal_domains( const libMesh::MeshBase & mesh )
{
  set<string> thermal;
  CFG.strset( "mesh", "thermal", thermal );
  dlog(1) << "LIDO SET - POROELASTIC::THERMAL:";
  for ( auto s : thermal ) dlog(1) <<"   " << s;
  mesh.print_info();
  set<uint> ret;
  for ( auto & s : thermal ) 
  {
    uint sid = mesh.get_id_by_name( s );
    if ( sid == libMesh::Elem::invalid_subdomain_id ) flog << "Dominio '" << s << "' nao existe na malha. Arquivo de configuraçao invalido."; 
    ret.insert(sid);
  }

  return ret;
}
/**
 *
 * Converte o nome dos dominios nos IDs dos dominios (FLUXO)
 *
 */
set<uint> MeshConfig::flow_domains( const libMesh::MeshBase & mesh )
{
  set<string> flow;
  CFG.strset( "mesh", "flow", flow );
  dlog(1) << "LIDO SET - POROELASTIC::FLOW:";
  for ( auto s : flow ) dlog(1) <<"   '" << s << "'";
  mesh.print_info();

  set<uint> ret;
  for ( auto & s : flow ) 
  {
    uint sid = mesh.get_id_by_name( s );
    if ( sid == libMesh::Elem::invalid_subdomain_id ) 
    {
      set<subdomain_id_type> sds; mesh.subdomain_ids(sds);
      dlog(1) << "Subdominios conhecidos:";
      for ( auto sd : sds ) dlog(1) << "    " << sd << ": " << mesh.subdomain_name(sd);
      flog << "Dominio '" << s << "' nao existe na malha. Arquivo de configuraçao invalido."; 
    }
    ret.insert(sid);
  }

  return ret;
}
