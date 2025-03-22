
#include "solver/BC.h"
#include "config/BCConfig.h"
#include "libmesh/elem.h"
#include "libmesh/boundary_info.h"

/**
 *
 * Initialize an empty boundary constriant object.
 *
 */
BC::BC( const System & sys_ ) : system(sys_), config( system.name() ), time(-999), reftime(-999)
{

}

/**
 * Updates the BC with information from the BCConfig at time t
 */
void BC::update( double t ) 
{
  time = t;

  /** Test if the reference time has changed for the new time **/
  double rt = config.get_reftime( time );
  if ( rt == reftime ) {
    dlog(1) << "Reftime did not change. Continuing...";
    return; // nothing has changed
  }
  reftime = rt;
  // Clean structures
  bcmap.clear();
  dirichlet.clear();

  // Perform each part of the update
  _validate();
  _update_dirichlet();
  _update_stot();
}

/**
 *
 */
void BC::_validate()
{
  /** Check if all the sides are present in the mesh. Show warnings otherwise. **/
  const BCConfig::TimeEntry & timeentry = config.entry_by_time[reftime];

  set<string> seen_bnames; 
  const MeshBase & mesh = system.get_mesh();
  const BoundaryInfo & bi = mesh.get_boundary_info();

  /** Iterate over all element sides **/
  for (const auto & elem : mesh.active_element_ptr_range()) 
  for (auto side : elem->side_index_range()) 
  {
    vector<boundary_id_type> vec;
    bi.boundary_ids(elem, side, vec );
    /** If the side has a boundary id ... **/
    for ( auto bid : vec ) 
    {
      /** Is the this bid currently constrained? **/
      string bname = bi.get_sideset_name( bid );
      if ( ! timeentry.dbl_bcs.count( bname ) ) continue;
      seen_bnames.insert( bname );
    }
  }

  /** Check if all defined boundaries in BC have at least one element side in the mesh **/
  set<string> bnames;
  config.all_bnames( bnames );
  for ( auto bname : bnames ) 
  if ( ! seen_bnames.count( bname ) )
    wlog << "No element has boundary '" << bname << "'. Check the mesh or the boundary conditions.";
}


/**
 * Update the Dirichlet BCs (bring from BCConfig)
 */
void BC::_update_dirichlet()
{
  const BCConfig::TimeEntry & timeentry = config.entry_by_time[reftime];

  // Refresh data structures - can be only the local (would save some minor time)
  const MeshBase & mesh = system.get_mesh();
  const BoundaryInfo & bi = mesh.get_boundary_info();

  for ( auto & [ bname, vec ] : timeentry.dbl_bcs )
  for ( auto & item : vec )
  {
    string vname = item.vname;

    /** Validation of the variable name (Dirichlet constraints are on
        existing variables of the system) **/
    if ( ! system.has_variable( vname ) ) continue; 
    uint vid = system.variable_number( vname );

    double val = item.value;
    int bid = bi.get_id_by_name( bname );

    DirichletItem di( bid, vid, val, vname, bname );
    dirichlet.push_back( di );
  }
}

/**
 * Update the Total Stress BCs (bring from BCConfig)
 */
void BC::_update_stot()
{
  const BCConfig::TimeEntry & timeentry = config.entry_by_time[reftime];

  // Refresh data structures - can be only the local (would save some minor time)
  const MeshBase & mesh = system.get_mesh();
  const BoundaryInfo & bi = mesh.get_boundary_info();

  for ( auto & [ bname, item ] : timeentry.stot_bcs )
  {
    string vname = item.vname;

    const vector<vector<double>> &v = item.value;

    RealTensor val( v[0][0], v[0][1], v[0][2],
                    v[1][0], v[1][1], v[1][2],
                    v[2][0], v[2][1], v[2][2]  );
    int bid = bi.get_id_by_name( bname );

    STotItem si( bid, val, bname );
    stot.push_back( si );
  }
}

/**
 * Update the double BCs (bring from BCConfig)
 */
//void BC::_update_dbl()
//{
//  const BCConfig::TimeEntry & timeentry = config.entry_by_time[reftime];

//  set<string> seen_bnames; // For validation purposes

//  // Refresh data structures - can be only the local (would save some minor time)
//  const MeshBase & mesh = system.get_mesh();
//  const BoundaryInfo & bi = mesh.get_boundary_info();

//  /** Iterate over all element sides **/
//  for (const auto & elem : mesh.active_element_ptr_range()) 
//  for (auto side : elem->side_index_range()) {
//    vector<boundary_id_type> vec;
//    bi.boundary_ids(elem, side, vec );
//    /** If the side has a boundary id ... **/
//    for ( auto bid : vec ) 
//    {
//      /** Is the this bid currently constrained? **/
//      string bname = bi.get_sideset_name( bid );
//      if ( ! timeentry.dbl_bcs.count( bname ) ) continue;
//    
//      const BCConfig::ItemDbl & item = timeentry.dbl_bcs.at(bname);
//      
//      /** Validation of the variable name **/
//      if ( ! system.has_variable( item.vname ) ) 
//        flog << "Cannot add variable '" << item.vname << "' as in DBL BCS. Unexisting in System '" << system.name() << "'.";

//      /** Add the resolved ids to the bcmap **/
//      uint vid = system.variable_number( item.vname );
//      ElemSide es( elem->id() , side );
//      vector<Item> v = bcmap[es];
//      v.push_back( Item(vid, item.value) );
//      bcmap[es] = v;

//      seen_bnames.insert( bname );
//    }
//  }

//  /** Check if all defined boundaries in BC have at least one element side in the mesh **/
//  set<string> bnames;
//  config.all_bnames( bnames );
//  for ( auto bname : bnames ) 
//  if ( ! seen_bnames.count( bname ) )
//    wlog << "No element has boundary '" << bname << "'. Check the mesh or the boundary conditions.";
//}

/**
 * Output stream operators
 */
ostream& operator<<(ostream& os, const BC & m)
{
  os << "Current boundary condition (BC):" << endl;
  os << "   DOUBLE BCS:" << endl;
  for ( auto & [es, item_vec] : m.bcmap )
  {
    os << "      " << es << " =>  " ;
    for ( auto & item : item_vec ) os << setw(20) << item;
    os << endl;
  }
  os << "   DIRICHLET BCS:" << endl;
  os << m.dirichlet;
  os << "   STOT BCS:" << endl;
  os << m.stot;
  return os;
}
ostream& operator<<(ostream& os, const BC::Item & m)
{
  os << "Item (vid:" << m.vid << ", val:" << m.val << ")";
  return os;
}
ostream& operator<<(ostream& os, const BC::ElemSide & m)
{
  os << "ElemSide (eid:" << m.eid << ", side:" << m.side << ")";
  return os;
}
ostream& operator<<(ostream& os, const BC::DirichletItem & m)
{
  os << "DirichletItem (bid:" << m.bid << "(" << m.bname << "), vid:" << m.vid << "(" << m.vname << "), val:" << m.val << ") " ;
  return os;
}
ostream& operator<<(ostream& os, const BC::STotItem & m)
{
  os << "STotItem (bid:" << m.bid << "(" << m.bname << "), val:" <<endl << m.val << ") " ;
  return os;
}
ostream& operator<<(ostream& os, const vector<BC::DirichletItem> & m)
{
  for ( auto di : m ) { os << setw(30) << di << endl; }
  return os;
}
ostream& operator<<(ostream& os, const vector<BC::STotItem> & m)
{
  for ( auto di : m ) { os << setw(30) << di << endl; }
  return os;
}

/**
 * Enable the object to index a map
 */
bool BC::ElemSide::operator<(const BC::ElemSide & other) const
{
  if (eid == other.eid) return side < other.side;
  return eid < other.eid;
}
