
#include "solver/viscoplastic/VPBC.h"
#include "config/ModelConfig.h"
#include "config/BCConfig.h"
#include "timeloop/Timestep.h"
#include "util/OutputOperators.h"

#include "libmesh/elem.h"
#include "libmesh/boundary_info.h"

namespace solver {
namespace viscoplastic {

using config::MODEL;

/**
 *
 * Initialize an empty boundary constriant object.
 *
 */
BC::BC( const System & sys_ ) : 
              system(sys_), 
              config( MODEL->boundary_config ),
              time(-999), reftime(-999),
              temporal_bcs( config )
{ }


/**
 * Cleanup the generated dynamic datastructures
 */
void BC::_cleanup()
{
  for ( auto item : stotitem_ptrs ) delete(item);
  stot.clear();
  stotitem_ptrs.clear();

  for ( auto [ k, item ] : penalty_ptrs ) delete(item);
  penalty.clear();
  penalty_ptrs.clear();
}

/**
 * Updates the BC with information from the BCConfig at time t
 *
 * Returns true if reftime has changed from previous update.
 */
bool BC::update( double t ) 
{
  SCOPELOG(1);

  time = t;

  /** Test if the reference time has changed for the new time **/
  double rt = config.get_reftime( time );
  if ( rt == reftime ) {
    dlog(1) << "Reftime did not change. Continuing...";
    return false; // nothing has changed
  }
  reftime = rt;

  // Clear structures
  dirichlet.clear();

  // Perform each part of the update
  _validate();
  _update_dirichlet();
  _update_scalar();
  _update_stot();
  _update_penalty();

  return true;
}

/**
 *
 */
bool BC::update( const Timestep & ts ) { return update( ts.time ) ; }

/**
 *
 */
void BC::_validate()
{
  /** Check if all the sides are present in the mesh. Show warnings otherwise. **/
//  const BCConfig::TimeEntry & timeentry = config.entry_by_time[reftime];

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
      seen_bnames.insert( bname );
    }
  }

  /** Check if all defined boundaries in BC have at least one element side in the mesh **/
  set<string> bnames;
  config.all_bnames( bnames );
  for ( auto bname : bnames ) 
  if ( ! seen_bnames.count( bname ) )
    wlog << "No element has boundary '" << bname << "'. Check the mesh or the boundary conditions.";

  using util::operator<<;
  dlog(1) << "Seen boundary names in the mesh: " << seen_bnames;
  dlog(1) << "Constrained boundaries in the current BC object: " << bnames;
}

/**
 * Update the Scalar BCs (bring from BCConfig)
 */
void BC::_update_scalar()
{
  const BCConfig::TimeEntry & timeentry = config.entry_by_time[reftime];

  // Refresh data structures - can be only the local (would save some minor time)
  const MeshBase & mesh = system.get_mesh();
  const BoundaryInfo & bi = mesh.get_boundary_info();

  for ( auto & [ bname, vec ] : timeentry.scalar_bcs )
  for ( auto & item : vec )
  {
    string vname = item.vname;

    /** Validation of the variable name (Dirichlet constraints are on
        existing variables of the system) **/
    if ( ! system.has_variable( vname ) ) continue; 
    uint vid = system.variable_number( vname );

    string scalar_name = item.value;

    // The number of the scalar variable
    if ( ! system.has_variable( scalar_name ) ) 
      flog << "Trying to attach to a scalar variable that does not exist ('" << scalar_name << "')! This should have been checked in ModelReader! Something is wrong!";

    uint svid = system.variable_number( scalar_name );

    int bid = bi.get_id_by_name( bname );

    ScalarItem si( bid, vid, svid, scalar_name, vname, bname );
    scalar.push_back( si );
  }
}

/**
 * Update the Dirichlet BCs (bring from BCConfig)
 */
void BC::_update_dirichlet()
{
  SCOPELOG(1);

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
    if ( ! system.has_variable( vname ) ) 
    {
      dlog(1) << "System '" << system.name() << "' HAS NO variable named '" << vname << "'. Skipping constrain.";
      continue; 
    }
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
  stot.clear();
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

    // Create and register item
    STotItem * sitem = new STotItem( bid, val, bname );
    stotitem_ptrs.insert(sitem);

    // Register every element and sid with this stress item
    for (const auto & elem : mesh.active_element_ptr_range()) 
    for (auto side : elem->side_index_range()) 
    {
      vector<boundary_id_type> vec;
      bi.boundary_ids(elem, side, vec );
      for ( auto _bid : vec ) 
      if ( bid == _bid )
      {
        ElemSide es(elem->id(), side);
        stot[es] = sitem;
        break;
      }
    }
  }
}

/**
 * Update the Penalty BCs (bring from BCConfig)
 */
void BC::_update_penalty()
{
  penalty.clear();
  const BCConfig::TimeEntry & timeentry = config.entry_by_time[reftime];

  // Refresh data structures - can be only the local (would save some minor time)
  const MeshBase & mesh = system.get_mesh();
  const BoundaryInfo & bi = mesh.get_boundary_info();

  for ( auto & [ bname, vec_item ] : timeentry.penalty_bcs )
  for ( auto & item : vec_item )
  {
    string vname = item.vname;
    const string & pen_id = item.value;

    if ( ! config.penalty.count( pen_id ) ) flog << "Cannot fine penalty entry '" << pen_id << "'. Should have been checked in BCConfig. Something is wrong.";
    const auto & pbc = config.penalty[pen_id];

    if ( ! system.has_variable( vname ) ) 
      flog << "Cannot find variable '" << vname << "' in the system '" << system.name() << "'.";
    uint vid = system.variable_number( vname );

    int bid = bi.get_id_by_name( bname );

    // Create and register item as needed
    PenaltyItem * pitem = 0;
    if ( penalty_ptrs.count(pen_id ) )
      pitem = penalty_ptrs[pen_id];
    else 
    {
      pitem = new PenaltyItem( bid, vid, pbc.K, pbc.value, pen_id );
      penalty_ptrs[pen_id] = pitem;
    }

    // Register every element and sid with this stress item
    for (const auto & elem : mesh.active_element_ptr_range()) 
    for (auto side : elem->side_index_range()) 
    {
      vector<boundary_id_type> vec;
      bi.boundary_ids(elem, side, vec );
      for ( auto _bid : vec ) 
      if ( bid == _bid )
      {
        ElemSide es(elem->id(), side);
        vector<PenaltyItem *> vp = penalty[es];
        vp.push_back( pitem );
        penalty[es] = vp;

        break;
      }
    }
  }
}

/**
 * Output stream operators
 */
ostream& operator<<(ostream& os, const BC & m)
{
  os << endl;
  os << "Current boundary condition (BC):" << endl;
  os << "   Time:" << m.time << setw(15) << "reftime:" << m.reftime << endl;
  os << "   DIRICHLET BCS:" << endl;
  os << m.dirichlet;
  os << "   SCALAR BCS:" << endl;
  os << m.scalar;
  os << "   PENALTY BCS:" << endl;
  os << m.penalty;
  os << "   STOT BCS:" << endl;
  os << m.stot;
  return os;
}
ostream& operator<<(ostream& os, const BC::ElemSide & m)
{ os << "ElemSide (eid:" << m.eid << ", side:" << m.side << ")"; return os; }
ostream& operator<<(ostream& os, const BC::DirichletItem & m)
{ os << "DirichletItem (bid:" << m.bid << "(" << m.bname << "), vid:" << m.vid << "(" << m.vname << "), val:" << m.val << ") " ; return os; }
ostream& operator<<(ostream& os, const BC::ScalarItem & m)
{ os << "ScalarItem (bid:" << m.bid << "(" << m.bname << "), vid:" << m.vid << "(" << m.vname << "), svid:" << m.svid << ", scalar_name:" << m.scalar_name ; return os; }
ostream& operator<<(ostream& os, const BC::PenaltyItem & m)
{ os << "PenaltyItem (bid:" << m.bid << "), vid:" << m.vid << ", pen_name: (" << m.pen_name << "), pen_K:" << m.pen_K << ", pen_val:" << m.pen_val ; return os; }
ostream& operator<<(ostream& os, const BC::STotItem & m)
{ os << "STotItem (bid:" << m.bid << "(" << m.bname << "), val:" <<endl << m.val << ") " ; return os; }
ostream& operator<<(ostream& os, const vector<BC::DirichletItem> & m)
{ for ( auto di : m ) { os << setw(30) << di << endl; } return os; }
ostream& operator<<(ostream& os, const vector<BC::ScalarItem> & m)
{ for ( auto si : m ) { os << setw(30) << si << endl; } return os; }
ostream& operator<<(ostream& os, const vector<BC::PenaltyItem *> & m)
{ for ( auto di : m ) { os << "(" << di << ") " << *di << endl; } return os; }

ostream& operator<<(ostream& os, const map<BC::ElemSide, BC::STotItem *> & m)
{
  for ( auto [ es, sitem ] : m )
  {
    os << setw(30) << es << " : ";
    os << "(" << sitem << ") " << *sitem << endl;
  }
  return os; 
}
ostream& operator<<(ostream& os, const map<BC::ElemSide, vector<BC::PenaltyItem *>> & m)
{ for ( auto [ es, vec_pitem ] : m ) os << setw(30) << es << " : " << vec_pitem; return os; }

/**
 * Enable the object to index a map
 */
bool BC::ElemSide::operator<(const BC::ElemSide & other) const
{
  if (eid == other.eid) return side < other.side;
  return eid < other.eid;
}

}} // ns
