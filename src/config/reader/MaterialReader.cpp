
#include "config/reader/MaterialReader.h"
#include "config/reader/ReaderRegex.h"
#include "config/MaterialConfig.h"

#include "util/File.h"
#include "util/String.h"
#include <iomanip>
#include <regex>
#include <set>

using namespace MRDEF;
using namespace harpy_string;
set<string> KNOWN_VAR_TYPES = { "CON", "FILE" };

/**
 *
 */
MaterialReader::MaterialReader( MaterialConfig & _config ) : config(_config) 
{
  check_files();
  parse_material_file();


}

/**
 *
 */
void MaterialReader::check_files()
{
  string model_dir = config.model_dir;
  string filename = config.filename;;
  if ( ! dir_exists( model_dir ) )
    flog << "Cannot find model directory '" << model_dir << "'!";

  if ( ! file_exists( filename ) )
    flog << "Cannot find material file '" << filename << "'!";

  dlog(1) << "Model directory and material file ok! We are good to go! (" << model_dir << ", " << filename << ")";
}

/**
 *  Reads the MATERIAL file. No need for a state machine as this is a linear file.
 */
void MaterialReader::parse_material_file()
{
  string filename = config.filename;;
  ifstream file(filename);

  if ( !file.is_open() ) flog << "Cannot open file " << filename << ". Should have been resolved in check_files. What's wrong?";

  ln = 0;

  smatch match;
  while (  getline(file, line)  ) 
  {
    ln++;
    line = remove_comments_and_trim( line );

    if ( next_state() ) continue;

    switch (current_state)
    {
      case State::INITIAL: { engine_state(); break; }
      case State::FEM: { fem_state(); break; }
      default: break;
    }
  }
}

/**
 *   Updates the machine state. Returns true if it was updated.
 *
 */
bool MaterialReader::next_state() 
{
  smatch match;

  CIMap<State> nextState = {
    { "engine",    State::ENGINE  },
    { "fem",       State::FEM     }
  };

  // Resets the state
  if (regex_match(line, RE_EMPTY)) { current_state = State::INITIAL; return true; }

  // New section. Changes the state
  if ( ! regex_search(line, match, RE_SEC_NAME) ) return false;

  string sec = match[1];

  if ( ! nextState.count(sec) ) flog << "Cannot find next state for section key '" << sec << "' (line " << ln << ").";

  // Save
  current_state  = nextState[ sec ];

  return true;
}

/**
 *
 */
void MaterialReader::engine_state() 
{
  smatch match;
  string vname;
  if ( regex_search( line, match, RE_STR_STR_STR ) ) 
    reg_param_str( match[1], match[2], match[3] );

  else if ( regex_search( line, match, RE_STR_NUM ) ) 
    reg_param_dbl( match[1], "CON", stod(match[2]) );

  else if ( regex_search( line, match, RE_STR_STR_NUM ) ) 
    reg_param_dbl( match[1], match[2], stod(match[3]) );

}

/**
 *
 */
void MaterialReader::fem_state() 
{
  smatch match;
  string vname;
  if ( regex_search( line, match, RE_STR_STR ) ) 
  {
    string key = match[1], val = match[2];
    if ( iequals( key, "type" ) ) config.fem_type = val;
    else
    if ( iequals( key, "family") ) config.fem_family = val;
    else
    if ( iequals( key, "order") ) config.fem_order = val;
    else
      flog << "Unknown key '" << key << "' in material parsing, FEM section.";
  }
}


/**
 *
 */
void MaterialReader::reg_param_str( string vname, string type, string val )
{
  if ( ! KNOWN_VAR_TYPES.count(type) ) flog << "Unkwown variable type at '" << filename << "' (" << ln << "): " << type;
  if ( iequals( type, "file" ) ) config.file_param( vname ) = val;
}

/**
 *
 */
void MaterialReader::reg_param_dbl( string vname, string type, double val )
{
  if ( iequals( type, "con" ) )  config.con_param(vname) = val ;
}
