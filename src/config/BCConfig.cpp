#include "config/BCConfig.h"
#include "config/Config.h"
#include "util/OutputOperators.h"
using namespace rapidjson;

/**
 *
 * Builds the datastructure from the json
 *
 */
BCConfig::BCConfig()  {}

/**
 *
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
  }
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

//{
////  build_scalars();
////  build_initial();
////  build_penalty();
////  build_bcs();
//}

///**
// *
// *
// */
//void BCConfig::build_scalars()
//{
//  SCOPELOG1(1);

//  if ( CFG.exists( sys_name, "scalars" ) )
//    for ( auto & V : CFG._rj[sys_name.c_str()]["scalars"].GetArray() )
//      scalars.insert(V.GetString());
//}

///*
// *
// * Returns the time imediatelly before time.
// *
// */
//double BCConfig::get_reftime( double time ) 
//{
//  double reftime = -999;
//  for ( const auto & [ t, e ] : entry_by_time )
//  {
//    if ( t > time ) break; 
//    reftime = t;
//  }
//  return reftime;
//}

///**
// *
// *
// *
// */
//void BCConfig::build_initial() 
//{
//  SCOPELOG1(1);

//  // Undefined initial conditions? No problem... Move on with empty arrays.
//  if ( ! CFG.exists( sys_name, "initial" ) ) 
//  {
//    dlog(1) << "Nao encontrei o namespace '"<<sys_name<<".initial' no json!" ;
//    return;
//  }

//  const Value& VV = CFG._rj[sys_name.c_str()]["initial"];

//  for (Value::ConstMemberIterator Vi = VV.MemberBegin(); Vi != VV.MemberEnd(); ++Vi) {
//    string var = Vi->name.GetString();
//    set<string> VARS = { "UX","UY","UZ","P","T" };

//    // O nome da variavel nao foi solicitado
//    if ( ! VARS.count(var) ) continue;

//    if (! Vi->value.IsNumber() ) flog << "Condicao inicial da variavel '"<<var<<"' nao eh numerica?";

//    double val = Vi->value.GetDouble();
//    initial_by_vname[var] = val;
//  }
//}


///**
// *
// * Constroi os penalties.
// *
// */
//void BCConfig::build_penalty() 
//{
//  SCOPELOG1(1);
//  if ( CFG.exists( sys_name.c_str(), "penalty" ) )
//  {
//    auto & X = CFG._rj[sys_name.c_str()]["penalty"];
//    for (Value::ConstMemberIterator Xi = X.MemberBegin(); Xi != X.MemberEnd(); ++Xi) {
//      string pname = Xi->name.GetString();
//      auto & Y = Xi->value;

//      PenaltyBC bc;

//      for (Value::ConstMemberIterator Yi = Y.MemberBegin(); Yi != Y.MemberEnd(); ++Yi) {
//        string var = Yi->name.GetString();
//        double val = Yi->value.GetDouble();

//        if ( var == "strength" ) bc.K = val;
//        if ( var == "value" )    bc.value = val;
//      }
//      penalty[pname] = bc;
//    }
//  }
//}

///**
// *
// *
// */
////void BC::scalar_bcs( BCMap<uint> & ret, set<uint> vars, const System & sys ) const
////{
////  for ( const auto & [bid,vec] : scalar_by_bid )
////    for ( const auto & [vid,rname] : vec )
////      if ( vars.count( vid ) )
////      {
////        if ( ! sys.has_variable( rname ) )
////        {
////          flog << "System does not have the scalar variable '"<< rname <<"'. Is it a penalty variable? Ignoring...";
////          continue;
////        }

////        ret.add( bid, vid, sys.variable_number(rname) );
////      }
////}

///**
// *
// *
// */
////void BC::penalty_bcs( BCMap<PenaltyBC> & ret, set<uint> vars  ) const
////{
////  for ( const auto & [bid,vec] : penalty_by_bid )
////  for ( const auto & [vid,pen] : vec )
////  if ( vars.count( vid ) )
////    ret.add( bid, vid, pen );
////}

///**
// *
// *
// */
////void BC::double_bcs( BCMap<double> & ret, set<uint> vars ) const
////{
////  for ( const auto & [bid,vec] : dbls_by_bid )
////    for ( const auto & [vid,val] : vec )
////      if ( vars.count( vid ) )
////        ret.add( bid, vid, val );
////}

/////**
//// *
//// *
//// */
////void BC::double_bcs( map<boundary_id_type, double> & ret, uint _vid ) const
////{
////  for ( const auto & [bid,vec] : dbls_by_bid )
////    for ( const auto & [vid,val] : vec )
////      if ( vid == _vid )
////        ret.insert( { bid, val } );
////}



///**
// *
// * Constroi as estruturas de dados
// *
// * estrutura:
// * CFG => poroelastic :{ 
// *             'drained' : bool,   // the default drained condition
// *             boundary_conditions : 
// *             {
// *                  'time' : dbl,
// *                   'drained' : bool, // drained condition of the timestep
// *                   'bc' : {
// *                       'bname':
// *                             {
// *                             "var1" : dbl,
// *                             "var2": dbl, ...
// *                             }
// *                    },
// *                    'flow' : {
// *                       'bname' : dbl
// *                    }
// *            }
// */
//void BCConfig::build_bcs() 
//{
//  SCOPELOG1(1);

//  // No boundary conditions? No problem. Move on with empty arrays
//  if ( ! CFG.exists( sys_name, "boundary_conditions" ) ) 
//  {
//    dlog(1) << "Namespace '"<<sys_name<<".boundary_conditions' not found in json." ;
//    return;
//  }
//  const Value& dbc_config = CFG._rj[sys_name.c_str()]["boundary_conditions"];

//  // Iterate BCs
//  for ( auto & V : dbc_config.GetArray() ) 
//  {
//    if ( ! V.HasMember("time") ) flog << "No json, o " << sys_name << "::boundary_conditions nao achei a chave 'time'.";
//    double t = V["time"].GetDouble();
//    TimeEntry entry = entry_by_time[t];

//    CFG.bln(sys_name, "drained", entry.drained, entry.drained );             // The global drained definition
//    if ( V.HasMember("drained") ) entry.drained = V["drained"].GetBool();      // override if exist in the BC

//    if ( V.HasMember("temperature") ) { entry.has_temperature = true; entry.temperature = V["temperature"].GetDouble(); }
//    if ( V.HasMember("pressure") ) { entry.has_pressure = true; entry.pressure = V["pressure"].GetDouble(); }

//    /* BCs */
//    if ( V.HasMember("bc") )
//    {
//      auto & X = V["bc"];
//      // X: { bname -> { v1:d, v2:d, ... } }
//      for (Value::ConstMemberIterator Xi = X.MemberBegin(); Xi != X.MemberEnd(); ++Xi) 
//      {
//        auto bname = Xi->name.GetString();
//        auto & Y = Xi->value;

//        // Y: { v1:d, v2:d, ... }
//        for (Value::ConstMemberIterator Yi = Y.MemberBegin(); Yi != Y.MemberEnd(); ++Yi) {
//          string var = Yi->name.GetString();

//          /* Tensors */
//          set<string> STOT_VARS = {"SXX","SYY","SZZ","SXY","SXZ","SYZ"};
//          if ( STOT_VARS.count(var) ) 
//          {
//            if (! Yi->value.IsNumber() ) flog << "Nao suporta scalar para valores tensoriais.";

//            // Fetch or create item or create
//            ItemTensor item = entry.stot_bcs[bname];
//            item.vname = "STOT";
//            item.bname = bname;

//            double val = Yi->value.GetDouble();
//            if ( var == "SXX" )       item.value[0][0] = val;
//            else if ( var == "SYY" )  item.value[1][1] = val;
//            else if ( var == "SZZ" )  item.value[2][2] = val;
//            else if ( var == "SXY" )  { item.value[0][1] = val; item.value[1][0] = val;   } 
//            else if ( var == "SXZ" )  { item.value[0][2] = val; item.value[2][0] = val;   }
//            else if ( var == "SYZ" )  { item.value[1][2] = val; item.value[2][1] = val;   }
//            else flog << "Unknown variable named '" << var << "'";

//            entry.stot_bcs[bname] = item; // update datastruct
//          }

//          /* Doubles, scalars and penalties */
//          set<string> DBL_VARS = {"UX","UY","UZ","P","T","Q"};
//          if ( DBL_VARS.count(var) ) 
//          {
//            if (Yi->value.IsNumber() ) 
//            {
//              ItemDbl item;
//              item.vname = var;
//              item.bname = bname;
//              item.value = Yi->value.GetDouble();
//              auto & vec = entry.dbl_bcs[bname];
//              vec.push_back(item);
//              entry.dbl_bcs[bname] = vec;
//            }

//            // SCALARS
//            if (Yi->value.IsString() ) 
//            {
//              string scalar_name = Yi->value.GetString();
//              ItemStr item;
//              item.vname = var;
//              item.bname = bname;
//              item.value = scalar_name;

//              if ( scalars.count( scalar_name ) )
//              {
//                auto & vec = entry.scalar_bcs[bname];
//                vec.push_back(item);
//                entry.scalar_bcs[bname] = vec;
//              } 
//              else if ( penalty.count( scalar_name ) )
//              {
//                auto & vec = entry.penalty_bcs[bname];
//                vec.push_back(item);
//                entry.penalty_bcs[bname] = vec;
//              }
//              else 
//                flog << "Cannot find '" << scalar_name << "' in scalar nor penalty variable list.";
//            }
//          }
//        }
//      }

//      // Register the entry and continue
//      entry_by_time[t] = entry;

//    } // bc
//  }
//}



/**
 *
 *
 *
 */
void BCConfig::all_bnames( set<string> & ret ) 
{
  for (const auto& [ts, entry] : entry_by_time ) 
  {
    for ( const auto & [bname, item] : entry.dbl_bcs ) ret.insert( bname );
    for ( const auto & [bname, item] : entry.scalar_bcs ) ret.insert( bname );
  }
}

/**
 *
 * Dumps the boundary conditions read
 *
 */
ostream& operator<<(ostream& os, const BCConfig & m)
{
  os << endl;
  os << "|-----------------------------------------------" << endl;
  os << "|   BOUNDARY CONDITIONS     " << endl;
  os << "|-----------------------------------------------" << endl;
  os << "| SCALARS: ";
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
ostream& operator<<(ostream& os, const BCConfig::TimeEntry & m) 
{
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

ostream& operator<<(ostream& os, const BCConfig::ItemTensor & m)
{
  const vector< vector<double> > & v = m.value;
  os << endl;
  os << "            [ [" << setw(10) << v[0][0] << "," << setw(10) << v[0][1] << "," << setw(10) << v[0][2] << "] ," << endl;
  os << "              [" << setw(10) << v[1][0] << "," << setw(10) << v[1][1] << "," << setw(10) << v[1][2] << "] ," << endl;
  os << "              [" << setw(10) << v[2][0] << "," << setw(10) << v[2][1] << "," << setw(10) << v[2][2] << "] ]" << endl;
  return os;
}

//string BC::dump( const System & sys ) const 
//{
//  const MeshBase & mesh = sys.get_mesh();
//  const BoundaryInfo & bi = mesh.get_boundary_info();

//  ostringstream os;
//  os << endl;
//  os << "| >>> ts= " << ts << "<<< drained?" << ( drained ? "TRUE" : "FALSE" ) << endl;
//  os << "| SIGTOT:" << endl << sigtot_by_bid.dump(sys) << endl;
//  os << "| DOUBLES:" << endl << dbls_by_bid.dump(sys) << endl;
//  os << "| HEAT:" << endl << heat_by_bid.dump(sys) << endl;
//  os << "| FLOW:" << endl; 
//  for (const auto& [bid, dbl] : flow_by_bid) {
//    string bname    = bi.get_sideset_name( bid );
//    os << "| BID: " << bname << "(" << bid << ") => " << dbl << endl;
//  }


//  return os.str();
//}

/**
 *
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
