
#include "config/ModelReader.h"
#include "config/ModelConfig.h"
#include "config/ReaderRegex.h"
#include "config/SystemConfig.h"

#include "util/File.h"
#include <iomanip>
#include <regex>
#include <set>


using namespace harpy_string;
using namespace MRDEF;

/**
 *  
 */
ModelReader::ModelReader( ModelConfig & config_ ) : config(config_) 
{
  check_files();
  parse_model_file();

  for ( auto & [ sysname, syscfg ] : config.system_cfgid )
    config.systems.emplace( 
                    sysname, 
                    SystemConfig(config.model_dir, sysname, syscfg) 
              );
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
    line = remove_comments_and_trim( line );

    if ( next_state() ) continue;

    // Machine. This is the last part of the loop
    switch (current_state) 
    {
      case State::TIMESTEP: { timestep_state(); break; }
      case State::SYSTEM: { system_state() ; break; }
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
    { "system", State::SYSTEM },
  };

  // Resets the state
  if (regex_match(line, emptyRE)) { current_state = State::INITIAL; return true; }

  dlog(1) << tok_op;

  // New section. Changes the state
  if ( ! regex_search(line, match, sectionRE) ) return false;
  string sec = match[1];
  if ( ! nextState.count(sec) ) flog << "Cannot find next state for section key '" << sec << "' (line " << ln << "'.";
  current_state = nextState[ sec ];
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
