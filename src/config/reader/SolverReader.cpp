
#include "config/reader/SolverReader.h"
#include "config/SolverConfig.h"
#include "config/reader/ReaderRegex.h"

#include "util/File.h"
#include <iomanip>
#include <regex>
#include <set>

namespace config {
namespace reader {

using namespace util;
using namespace MRDEF;

/**
 *  
 */
SolverReader::SolverReader( SolverConfig & config_ ) : config(config_) 
{
  check_files();
  parse_sys_file();

//  dlog(1) << *this;

//  Set the right configuration
  string curr_cfg = config.sys_cfg;

  if ( iequals( curr_cfg , "external" ) ) 
  {
    config.external_file.check(); // This shall fail if not completely defined
    return;
  }

  if ( ! all_mat_cfgs.count( curr_cfg ) )
    flog << "System configuration '" << curr_cfg << "' not found for system '" << config.sys_name << "'.";

  // Select the material config
  config.mat_config_by_name = all_mat_cfgs[ curr_cfg ];
}

/**
 *
 */
string SolverReader::abs_filepath( string rel_filepath, bool check_exists )
{
  string ret = config.model_dir + "/" + rel_filepath;

  if ( check_exists )
  if ( ! file_exists( ret ) )
    flog << "Cannot find file '" << ret << "'!";

  return ret;
}

/**
 *
 */
void SolverReader::check_files()
{
  string model_dir = config.model_dir;
  string sys_file = config.sys_file;
  
  if ( ! dir_exists( model_dir ) )
    flog << "Cannot find model directory '" << model_dir << "'!";

  if ( ! file_exists( sys_file ) )
    flog << "Cannot find model file '" << sys_file << "'!";

  dlog(1) << "Model directory and system file ok! We are good to go! (" << model_dir << ", " << sys_file << ")";
}

/**
 *  Reads the MODEL file: this is the main STATE MACHINE for the MODEL file.
 *
 *  This function holds the state machine definitions.
 *  Regex and other stuff, add in the MRDEF namespace for clarity.
 *
 */
void SolverReader::parse_sys_file()
{
  ifstream file(config.sys_file);
  if ( !file.is_open() ) flog << "Cannot open file " << config.sys_file << ". Should have been resolved in check_files. What's wrong?";

  current_state = State::INITIAL;
  ln = 0;

  while (  getline(file, line)  ) 
  {
    ln++;
    line = remove_comments_and_trim( line );

    if ( next_state() ) continue;

    // Machine. This is the last part of the loop
    switch (current_state) 
    {
      case State::EXTERNAL:  {   external_state();      break;    }
      case State::CONFIG:    {   config_state();      break;    }
      case State::NUMERICAL: {   numerical_state();   break;    }
      case State::FEM:       {   fem_state();         break;    }
      default: break;
    } 
  } // machine loop
}

/**
 *   Updates the machine state. Returns true if it was updated.
 *
 */
bool SolverReader::next_state() 
{
  smatch match;

  CIMap<State> nextState = {
    { "external", State::EXTERNAL },
    { "config", State::CONFIG },
    { "mesh",   State::MESH },
    { "fem",   State::FEM },
    { "numerical",   State::NUMERICAL }
  };

  set<State> namedStates = { State::CONFIG, State::MESH }; // These require a name following the section

  // Resets the state
  if (regex_match(line, RE_EMPTY)) { current_state = State::INITIAL; return true; }

  // New section. Changes the state
  if ( ! regex_search(line, match, RE_SEC_NAMEOPT) ) return false;

  string sec = match[1];
  if ( ! nextState.count(sec) ) flog << "Cannot find next state for section key '" << sec << "' (line " << ln << ").";
  current_state  = nextState[ sec ];

  // New section. Changes the state
  if ( ! regex_search(line, match, RE_SEC_NAME) ) 
  {
    if ( ! namedStates.count( current_state ) ) return true; // all good
    else flog << "Keyword '" << sec << "' requires a name at solver definition (line " << ln << ").";
  }

  // Save
  string name = match[2];

  // Register as needed
  if ( current_state == State::CONFIG ) curr_sys_cfg = name;

  if ( current_state == State::MESH ) 
  {
    // Single line state.
    config.mesh_filename = config.model_dir + "/" + name;
    current_state = State::INITIAL;
    dlog(1) << "Read mesh filename: '" << config.mesh_filename << "'.";
  }

  return true;
}


/**
 *
 *
 */
void SolverReader::set_grid_type( string str ) { 
  using enum GRID_TYPE;

  if ( iequals( str, "radial" ) ) 
    config.external_file.grid_type = RADIAL;

  else flog << "Unknown grid type "<< str <<".";
}

/**
 * Parse a line in the state
 */
void SolverReader::external_state()
{
  dlog(1) << "Processing EXTERNAL FILE state line '" << line << "' ";
  smatch match;

  //
  if ( regex_search( line, match, RE_STR_STR ) ) 
  {

    string key = match[1];
    string val = match[2];

    if ( iequals( key, "file" ) )       { config.external_file.filename = abs_filepath(val); return; }
    if ( iequals( key, "grid" ) )       { set_grid_type(val);                                return; }
    if ( iequals( key, "min_radius" ) ) { config.external_file.min_radius = stod(val);       return; }
  }

  //
  if ( regex_search( line, match, RE_STR_POINT ) ) 
  {

    string key = match[1];
    double x = stod( match[2] ) , 
           y = stod( match[3] ) ,
           z = stod( match[4] );

    if ( iequals( key, "origin" ) ) { config.external_file.grid_origin = libMesh::Point( x, y, z ); return; }
      flog << ">>" << key << "<<";
  }

  flog << "Unrecognized format at " << config.sys_file << "(" << ln << "), EXTERNAL section. Line: " << line;

}


/**
 * Parse a line in the state
 */
void SolverReader::config_state()
{
  dlog(1) << "Processing CONFIG state line '" << line << "' ";
  smatch match;

  if ( regex_search( line, match, RE_STR_STR_STROPT ) ) 
  {
    string subdom = match[1];
    string material = match[2];
    string conf = match[3];
    if ( ! conf.length() ) conf = curr_sys_cfg;

    all_mat_cfgs[curr_sys_cfg].emplace( 
                   subdom,
                   SolverConfig::MatNameAndCfg( material, conf ) );
    return;
  }

  flog << "Unrecognized format at " << config.sys_file << "(" << ln << "), CONFIG section. Line: " << line;

}

/**
 * Parse a line in the state
 */
void SolverReader::numerical_state()
{
  dlog(1) << "Processing numerical state line '" << line << "' ";
  auto & num = config.numerical;

  smatch match;

  // We have a numerical parameter?
  if ( regex_search( line, match, RE_STR_NUM ) ) 
  {
    string vname = match[1];
    double val = stod( match[2] ); 
    to_lower(vname);
    
    if ( vname == "ls_atol" ) num.ls_atol = val;
    if ( vname == "ls_rtol" ) num.ls_rtol = val;

    return;
  } 
  // We have a string parameter?
  else if ( regex_search( line, match, RE_STR_STR ) ) 
  {
    string vname = match[1];
    string val = match[2];
    to_lower(vname);

    if ( vname == "ls_pc" )    num.ls_pc = val;
    if ( vname == "ls_rtol" )  num.ls_solver = val;

    return;
  } 

  flog << "Unrecognized format at " << config.sys_file << "(" << ln << "), MATERIAL section. Line: " << line;

}

/**
 *
 */
void SolverReader::fem_state() 
{
  SCOPELOG(1);

  using util::to_upper;

  const set<string> KNOWN_FE_TYPE = { "CONTINUOUS", "DISCONTINUOUS" };
  const set<string> KNOWN_FE_FAMILY = { "LAGRANGE", "L2_LAGRANGE", "ENRICHED_GALERKIN" };

  smatch match;
  string vname;
  if ( regex_search( line, match, RE_STR_STR_STR ) ) 
  {
    string key = match[1], var = match[2], val = match[3];
    to_upper(val); // upper case

    using FEMSpec = SolverConfig::FEMSpec;
    FEMSpec & fem = config.fem_by_var[var];
    //
    if ( iequals( key, "type" ) )
    {
      if (! KNOWN_FE_TYPE.count(val) ) flog << "Unknown FE type '" << val << "'.";
      fem.type = val;
    }
    //
    else if ( iequals( key, "family") )   
    {
      if (! KNOWN_FE_FAMILY.count(val) ) flog << "Unknown FE Family '" << val << "'.";
      fem.family = val;
    }
    //
    else if ( iequals( key, "order") )     fem.order = val;
    //
    else if ( iequals( key, "implicit") )  {
      if ( ! regex_search( val, match, RE_NUM ) ) flog << "Invalid value for IMPLICIT. Must be a number.";
      fem.implicit = stod(val);
    }
    else flog << "Unknown key '" << key << "' in material parsing, FEM section.";
  }
}


/**
 *
 *
 */
ostream& operator<<(ostream& os, const SolverReader & m)
{
  os << endl;
  os << "SolverReader for '" << m.config.sys_name  << "':" << endl;
  os << "   All material configurations" << endl;
  os << endl;
  for ( auto & [ sys_cfg , vm ] : m.all_mat_cfgs ) 
  {
    os << "------------------------------------------------------------------------------------------" << endl;
    os << "      Configuration'" << sys_cfg << "':" << endl << vm;
  }
  return os;
}


}} // ns
