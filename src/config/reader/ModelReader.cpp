
#include "config/reader/ModelReader.h"
#include "config/reader/ReaderRegex.h"
#include "config/ModelConfig.h"
#include "config/SolverConfig.h"

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
ModelReader::ModelReader( ModelConfig & config_ ) : config(config_), current_time(-1)
{
  check_files();
  parse_model_file();

  // Reads all the solvers
  for ( auto & [ sysname, syscfg ] : config.system_cfgid )
    config.solvers.emplace( 
                    sysname, 
                    SolverConfig(config.model_dir, sysname, syscfg) 
              );
  
  // Reads all the materials
  set<SolverConfig::MatNameAndCfg> mats;
  for ( auto & [ sysname, syscfg ] : config.solvers )
  for ( auto & [ matname, matcfg ] : syscfg.mat_config_by_name )
    mats.insert( matcfg );
  for ( auto & mat : mats )
    config.materials.emplace( MaterialConfig( config.model_dir, mat.name, mat.cfg ) );
}


/**
 *
 */
void ModelReader::check_files()
{
  string model_dir = config.model_dir;
  string model_file = config.model_file;
  
  if ( ! dir_exists( model_dir ) )
    flog << "Cannot find model directory '" << model_dir << "'!";

  if ( ! file_exists( model_file ) )
    flog << "Cannot find model file '" << model_file << "'!";

  dlog(1) << "Model directory and file ok! We are good to go! (" << model_dir << ", " << model_file << ")";
}

/**
 *  Reads the MODEL file: this is the main STATE MACHINE for the MODEL file.
 *
 *  This function holds the state machine definitions.
 *  Regex and other stuff, add in the MRDEF namespace for clarity.
 *
 */
void ModelReader::parse_model_file()
{
  ifstream file(config.model_file);
  if ( !file.is_open() ) flog << "Cannot open file " << config.model_file << ". Should have been resolved in check_files. What's wrong?";

  current_state = State::INITIAL;
  ln = 0;

  while (  getline(file, line)  ) 
  {
    ln++;
    trim_right(line);
    line = remove_comments_and_trim( line );

    if ( next_state() ) continue;

    // Machine. This is the last part of the loop
    switch (current_state) 
    {
      case State::TIMESTEP:   { timestep_state(); break; }
      case State::SYSTEMLOOP: { system_state() ; break; }
      case State::INITIAL:    { time_state() ; break; }
      case State::EXTERNAL:   { external_state() ; break; }
      case State::TIME:       { time_state() ; break; }
      case State::SCALAR:     { scalar_state() ; break; }
      case State::PENALTY:    { penalty_state() ; break; }
      default: break;
    } 

     
  } // machine loop
}

/**
 *   Updates the machine state. Returns true if it was updated.
 */
bool ModelReader::next_state() 
{
  smatch match;

  CIMap<State> nextState = {
    { "timestep", State::TIMESTEP },
    { "systemloop", State::SYSTEMLOOP },
    { "initial", State::INITIAL },
    { "external", State::EXTERNAL },
    { "scalar", State::SCALAR },
    { "penalty", State::PENALTY },
    { "time", State::TIME },
  };

  // Resets the state
  if (regex_match(line, RE_EMPTY)) { current_state = State::INITIAL; return true; }

  // New section. Changes the state
  if ( ! regex_search(line, match, RE_SEC) ) return false;
  string sec = match[1];
  if ( ! nextState.count(sec) ) flog << "Cannot find next state for section key '" << sec << "' (line " << ln << "'.";
  current_state = nextState[ sec ];


  /* Capture id of the section */
  switch (current_state) 
  {
    /**/
    case State::SYSTEMLOOP: { config.systemloop    = match[2];  break; }
    /**/
    case State::TIME: 
      { 
        string timestr = match[2];
        if ( ! regex_match( timestr, RE_NUM ) )
          flog<<"Time must be a number (line " << ln << "). Read '" << match[2] << "'.";
        current_time = stod( timestr );
        break; 
      }
    /**/
    default: break;
  } 

  return true;
}

/**
 * Parse a line in the state
 */
void ModelReader::timestep_state()
{
  smatch match;

  // Known keys ( case insensitive )
  CISet timestepKeys = { "t0", "tmax", "dt0", "dtk", "max_steps" };

  if ( ! regex_search( line, match, RE_STR_NUM ) ) 
    flog << "Unrecognized format at MODEL(" << ln << "), TIMESTEP section. Line: " << line;

  string key = match[1], val = match[2];

  config.timestep.set( key, stod( val ) );
}
/**
 * Parse a line in the state
 */
void ModelReader::system_state()
{
  smatch match;

  if (regex_search(line, match, RE_STR_STR)) 
  {
    config.system_cfgid[ match[1] ] = match[2];
  }
}
/**
 * Parse a line in the state
 */
void ModelReader::time_state()
{
  smatch match;

  // BOUNDARY/DOMAIN    BOUNDARY_NAME    VAR_NAME   VALUE

  // scalar, penalty or numerical
  if ( regex_search( line, match, RE_STR_STR_STR_STR ) ) 
  {
    string type = match[1];
    string scope = match[2];
    string vname = match[3];
    string strval = match[4];

    // Validation
    CISet types = { "domain", "boundary" };
    if ( ! types.count( type ) ) flog << "Unknown constraint type '" << type << "'. MODEL, line " << ln << ".";

    BCConfig & bcconfig = config.boundary_config;

    // Numerical
    if ( regex_match( strval, RE_NUM ) ) 
    {
      double val = stod( strval );
      auto & te = bcconfig.entry_by_time[ current_time ];
      if ( iequals( type, "domain" ) )
        te.add_domain_bc( scope, vname, val );
      else if ( iequals( type, "boundary" ) )
        te.add_numerical_bc( scope, vname, val );
    } 
    // Penalty/scalar
    else 
    {
      if ( ! iequals( type, "boundary" ) ) 
        flog << "Scalars and penalty constrains are only supported in BOUNDARY constraints (not in DOMAIN ones!). MODEL, line " << ln << ".";

      auto & te = bcconfig.entry_by_time[ current_time ];
      if ( bcconfig.penalty.count( strval ) )
        te.add_penalty_bc( scope, vname, strval );
      else if ( bcconfig.has_scalar( strval ) )
        te.add_scalar_bc( scope, vname, strval );
      else 
        flog << "Cannot find '" << strval << "' in scalar nor penalty variable list (MODEL file, line " << ln << ").";
    }
  } else
    flog << "Unrecognized format at MODEL(" << ln << "), INITIAL section. Line: " << line;
}

/**
 * Parse a line in the state
 */
void ModelReader::external_state()
{
  smatch match;


  // Invalid definition!
  if ( ! regex_search( line, match, RE_STR_STR ) ) flog << "Unrecognized format at MODEL(" << ln << "), EXTERNAL section. Line: " << line;

  string vname    = match[1],
  filename = match[2];

  // This is the only supported variable for now.
  if ( ! iequals( vname, "sig0" ) ) flog << "Unsupported variable at MODEL (" << ln << ").";

  // Check if file exists
  string full_fn = config.model_dir + "/" + filename;

  if ( ! file_exists( full_fn ) ) flog << "File '" << full_fn << "' does not exist. Review MODEL (" << ln << ") file.";

  config.sig0_file = full_fn;

  return; 
}



/**
 * Parse a line in the state
 */
void ModelReader::scalar_state()
{
  smatch match;

  BCConfig & bcconfig = config.boundary_config;

  // Invalid definition!
  if ( regex_search( line, match, RE_STR_STR ) ) flog << "Unrecognized format at MODEL(" << ln << "), SCALAR section. Line: " << line;

  // Full definition ( NAME FAMILY ORDER )
  if ( regex_search( line, match, RE_STR_STR_STR ) ) 
  { 
    string name = match[1], family=match[2], order=match[3];
    // Validation
    CISet CHECK_FMLY = {"LAGRANGE"};
    if (! CHECK_FMLY.count(family) ) flog << "Invalid FAMILY definition at MODEL(" << ln << "), SCALAR SECTION. Read: '" << family << "'.";
    CISet CHECK_ORDER = {"FIRST", "SECOND"};
    if (! CHECK_ORDER.count(order) ) flog << "Invalid ORDER definition at MODEL(" << ln << "), SCALAR SECTION. Read: '" << order << "'.";

    // Add.
    bcconfig.scalars.emplace( name, family, order ); 
    return; 
  }

  // Short definition ( NAME ) -- FAMILY and ORDER are defaults -- see ScalarVar class
  if ( regex_search( line, match, RE_STR ) ) { bcconfig.scalars.emplace( match[1] ); return; }

  flog << "Unrecognized format at MODEL(" << ln << "), SCALAR section. Line: " << line;
}
/**
 * Parse a line in the state
 */
void ModelReader::penalty_state()
{
  smatch match;
  if ( ! regex_search( line, match, RE_STR_NUM_NUM ) ) 
    flog << "Unrecognized format at MODEL(" << ln << "), PENALTY section. Line: " << line;

  string vname = match[1];
  double k = stod( match[2] );
  double val = stod( match[3] );

  BCConfig & bcconfig = config.boundary_config;
  bcconfig.penalty.emplace( vname, BCConfig::PenaltyBC( k, val ) );
}


}} // ns
