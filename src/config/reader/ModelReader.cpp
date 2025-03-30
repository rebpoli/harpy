
#include "config/reader/ModelReader.h"
#include "config/reader/ReaderRegex.h"
#include "config/ModelConfig.h"
#include "config/SolverConfig.h"

#include "util/File.h"
#include <iomanip>
#include <regex>
#include <set>


using namespace harpy_string;
using namespace MRDEF;

/**
 *  
 */
ModelReader::ModelReader( ModelConfig & config_ ) : config(config_), current_time(-1)
{
  check_files();
  parse_model_file();

  // Reads all the systems
  for ( auto & [ sysname, syscfg ] : config.system_cfgid )
    config.systems.emplace( 
                    sysname, 
                    SolverConfig(config.model_dir, sysname, syscfg) 
              );
  
  // Reads all the materials
  set<SolverConfig::MatConfig> mats;
  for ( auto & [ sysname, syscfg ] : config.systems )
  for ( auto & [ matname, matcfg ] : syscfg.material_config )
    mats.insert( matcfg );
  for ( auto & mat : mats )
    config.materials.emplace( mat.name, MaterialConfig( config.model_dir, mat.name, mat.cfg ) );
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
    dlog(1) << "LINE(" << ln << ") >>" << line;
    line = remove_comments_and_trim( line );

    if ( next_state() ) continue;

    // Machine. This is the last part of the loop
    switch (current_state) 
    {
      case State::TIMESTEP:   { timestep_state(); break; }
      case State::SYSTEMLOOP: { system_state() ; break; }
      case State::INITIAL:    { time_state() ; break; }
      case State::TIME:       { time_state() ; break; }
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
    { "scalar", State::SCALAR },
    { "penalty", State::PENALTY },
    { "time", State::TIME },
  };

  // Resets the state
  if (regex_match(line, RE_EMPTY)) { current_state = State::INITIAL; return true; }

  dlog(1) << tok_op;

  // New section. Changes the state
  if ( ! regex_search(line, match, sectionRE) ) return false;
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

  if ( ! timestepKeys.count(key) )
    flog << "Unrecognized key at MODEL("<< ln << "), TIMESTEP section. Key: " << key;

  config.timestep[key] = stod( val );
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
    string bname = match[2];
    string vname = match[3];
    string strval = match[4];
    dlog(1) << "   time_state: type(" << type << ") bname(" << bname << ") vname(" << vname << ") strval(" << strval << ")" ;

    BCConfig & bcconfig = config.boundary_config;

    // Numerical
    if ( regex_match( strval, RE_NUM ) ) 
    {
      double val = stod( strval );
      auto & te = bcconfig.entry_by_time[ current_time ];
      te.add_numerical_bc( bname, vname, val );
    } 
    // Penalty/scalar
    else 
    {
      auto & te = bcconfig.entry_by_time[ current_time ];
      if ( bcconfig.penalty.count( strval ) )
        te.add_penalty_bc( bname, vname, strval );
      else if ( bcconfig.scalars.count( strval ) )
        te.add_scalar_bc( bname, vname, strval );
      else 
        flog << "Cannot find '" << strval << "' in scalar nor penalty variable list (MODEL file, line " << ln << ").";
    }
  } else
    flog << "Unrecognized format at MODEL(" << ln << "), INITIAL section. Line: " << line;


}
/**
 * Parse a line in the state
 */
void ModelReader::scalar_state()
{
  smatch match;
  if ( ! regex_search( line, match, RE_STR ) ) 
    flog << "Unrecognized format at MODEL(" << ln << "), SCALAR section. Line: " << line;

  BCConfig & bcconfig = config.boundary_config;
  bcconfig.scalars.insert( match[1] );
}
/**
 * Parse a line in the state
 */
void ModelReader::penalty_state()
{
  smatch match;
  if ( regex_search( line, match, RE_STR_NUM_NUM ) ) 
    flog << "Unrecognized format at MODEL(" << ln << "), PENALRY section. Line: " << line;

  string vname = match[1];
  double k = stod( match[2] );
  double val = stod( match[3] );

  BCConfig & bcconfig = config.boundary_config;
  bcconfig.penalty.emplace( vname, BCConfig::PenaltyBC( k, val ) );
}
