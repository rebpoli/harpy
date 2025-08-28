#include "config/BCConfig.h"
#include "config/Config.h"
#include "util/OutputOperators.h"
#include "util/String.h"

namespace config {

/**
 *
 * Builds the datastructure from the  model.
 *
 */
BCConfig::BCConfig()  {}

/**
 *   A domain BC defines forces a variable to be constant in
 *   the timestep.
 */
void BCConfig::TimeEntry::add_domain_bc( string subdomain, string vname, double val )
{
  set<string> VARS = {"UX","UY","UZ","P","T"};
  if ( ! VARS.count(vname) ) flog << "Unknown variable for a domain BC: '" << vname << "'.";

  DomainBC dbc;
  dbc.vname = vname;
  dbc.subdomain = subdomain;
  dbc.value = val;

  auto & vec = domain_bcs[subdomain];
  vec.push_back(dbc);
  domain_bcs[subdomain] = vec;
}

/*
 *
 * Returns the time imediatelly before time.
 *
 */
double BCConfig::get_reftime( double time )
{
  double reftime = -999;
  for ( const auto & [ t, e ] : entry_by_time )
  {
    if ( t > time ) break;
    reftime = t;
  }
  return reftime;
}

/**
 *  Adds a numerical BC. Can be a stress, a displacement, pressure etc.
 */
void BCConfig::TimeEntry::add_numerical_bc( string bname, string vname, double val )
{
  // If this is a stress var ...
  set<string> STOT_VARS = {"SXX","SYY","SZZ","SXY","SXZ","SYZ"};
  if ( STOT_VARS.count(vname) ) 
  {
    ItemTensor item = stot_bcs[bname];
    item.vname = "STOT";
    item.bname = bname;
    if ( vname == "SXX" )       item.value[0][0] = val;
    else if ( vname == "SYY" )  item.value[1][1] = val;
    else if ( vname == "SZZ" )  item.value[2][2] = val;
    else if ( vname == "SXY" )  { item.value[0][1] = val; item.value[1][0] = val;   } 
    else if ( vname == "SXZ" )  { item.value[0][2] = val; item.value[2][0] = val;   }
    else if ( vname == "SYZ" )  { item.value[1][2] = val; item.value[2][1] = val;   }
    else flog << "Unknown variable named '" << vname << "'";
    stot_bcs[bname] = item; // update datastruct
    return;
  }

  // This is a Dirichlet
  set<string> DBL_VARS = {"UX","UY","UZ","P","T"};
  if ( DBL_VARS.count(vname) ) 
  {
    ItemDbl item;
    item.vname = vname;
    item.bname = bname;
    item.value = val;

    auto & vec = dbl_bcs[bname];
    vec.push_back(item);
    dbl_bcs[bname] = vec;
    return;
  }

  using util::operator<<;
  elog << "Known variables are: " << STOT_VARS << " and " << DBL_VARS;
  flog << "Unknown variable named '" << vname << "'";
}

/**
 *
 */
void BCConfig::TimeEntry::add_scalar_bc( string bname, string vname, string scalar_name )
{
  ItemStr item;
  item.vname = vname;
  item.bname = bname;
  item.value = scalar_name;
  auto & vec = scalar_bcs[bname];
  vec.push_back(item);
  scalar_bcs[bname] = vec;
}

/**
 *
 */
void BCConfig::TimeEntry::add_penalty_bc( string bname, string vname, string penalty_name )
{
  ItemStr item;
  item.vname = vname;
  item.bname = bname;
  item.value = penalty_name;
  auto & vec = penalty_bcs[bname];
  vec.push_back(item);
  penalty_bcs[bname] = vec;
}



/**
 *
 */
void BCConfig::all_bnames( set<string> & ret ) const
{
  for (const auto& [ts, entry] : entry_by_time )
  {
    for ( const auto & [bname, item] : entry.dbl_bcs ) ret.insert( bname );
    for ( const auto & [bname, item] : entry.scalar_bcs ) ret.insert( bname );
  }
}
/**
 *    Return true if a scalar named $name exists.
 */
bool BCConfig::has_scalar( string name ) const { return scalars.count( ScalarVar(name) ) > 0; }

/**
 *    Makes the ScalarVar object an indexer for maps and sets.
 *
 *    The object is only indexed by the _name_ variable, case insensitive.
 */
bool BCConfig::ScalarVar::operator<(const BCConfig::ScalarVar & other) const
{ util::ci_cmp c; return c( name, other.name ); }

/**
 *
 * Dumps the boundary conditions read
 *
 */
ostream& operator<<(ostream& os, const BCConfig & m)
{
  using util::operator<<;
  os << endl;
  os << "|-----------------------------------------------" << endl;
  os << "|   BOUNDARY CONDITIONS     " << endl;
  os << "|-----------------------------------------------" << endl;
  os << "| SCALARS: " << endl;
  os << m.scalars << endl;
  os << "|----------------------------------------------" << endl;
  os << "| INITIAL_BY_VNAME: " << endl;
  os << m.initial_by_vname << endl;
  os << "|----------------------------------------------" << endl;
  os << "| PENALTY: " << endl;
  os << m.penalty << endl;
  os << "|----------------------------------------------" << endl;
  os << "| ENTRY_BY_TIME: " << endl;
  os << m.entry_by_time << endl;
  os << "|----------------------------------------------" << endl;
  return os;
}

/**
 *
 */
ostream& operator<<(ostream& os, const map<double,BCConfig::TimeEntry> & m)
{
  for (const auto& [ts, item] : m ) {
    os << "|---------------------------" << endl;
    os << "| Time=" << ts << endl;
    os << "|---------------------------" << endl;
    os << item << endl;
    os << "|---------------------------" << endl;
  }
  return os;
}

/**
 *
 */
ostream& operator<<(ostream& os, const BCConfig::TimeEntry & m)
{
  os << "DOMAIN BCS:" << endl;
  for ( const auto & [subd, vec] : m.domain_bcs )
  for ( const auto & dbc : vec )
    os << "    " << subd << ": '" << dbc.vname << "'=" << dbc.value << endl;
  os << "DBL BCS:" << endl;
  for ( const auto & [bname, vec] : m.dbl_bcs )
  for ( const auto & item : vec )
    os << "    " << bname << ": '" << item.vname << "'=" << item.value << endl;
  os << "STR BCS:" << endl;
  for ( const auto & [bname, vec] : m.scalar_bcs )
  for ( const auto & item : vec )
    os << "    " << bname << ": '" << item.vname << "'=" << item.value << endl;
  os << "STOT BCS:" << endl;
  for ( const auto & [bname, item] : m.stot_bcs )
    os << "    " << bname << ": '" << item.vname << "'=" << item << endl;

  return os;
}

/**
 *
 */
ostream& operator<<(ostream& os, const BCConfig::ItemTensor & m)
{
  const vector< vector<double> > & v = m.value;
  os << endl;
  os << "            [ [" << setw(10) << v[0][0] << "," << setw(10) << v[0][1] << "," << setw(10) << v[0][2] << "] ," << endl;
  os << "              [" << setw(10) << v[1][0] << "," << setw(10) << v[1][1] << "," << setw(10) << v[1][2] << "] ," << endl;
  os << "              [" << setw(10) << v[2][0] << "," << setw(10) << v[2][1] << "," << setw(10) << v[2][2] << "] ]" << endl;
  return os;
}

/**
 *
 */
ostream& operator<<(ostream& os, const BCConfig::PenaltyBC & m)
{
  os << ": { 'K': " << m.K << ", 'value': " << m.value << " } (PenaltyBC)"; 
  return os;
}
ostream& operator<<(ostream& os, const map<string,BCConfig::PenaltyBC> & m)
{
  os << "All penalties (Map): {" << endl; 
  for ( const auto & [k,v] : m ) os << "    '" << k << "' : " << v << endl;
  os << "}" << endl;
  return os;
}

/**
 *
 */
ostream& operator<<(ostream& os, const set<BCConfig::ScalarVar> & m)
{
  os << "       ";
  os << setw(20) << "Name";
  os << setw(20) << "Family";
  os << setw(20) << "Order";
  os << endl;
  for ( auto & v : m ) os << v << endl;
  return os;
}
ostream& operator<<(ostream& os, const BCConfig::ScalarVar & m)
{
  os << "       ";
  os << setw(20) << m.name;
  os << setw(20) << m.family;
  os << setw(20) << m.order;
  return os;
}

}
