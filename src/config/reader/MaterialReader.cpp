
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
 *  Reads the MATERIAL file.
 */
void MaterialReader::parse_material_file()
{
  SCOPELOG(1);

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
    dlog(1) << "Parsing state: " << current_state;

    switch (current_state)
    {
      case State::POROTHERMOELASTIC: { porothermoelastic_state(); break; }
      case State::CREEP_CARTER: { creep_carter_state(); break; }
      default: { wlog << "Ignoring line " << ln << " >> " << line; break; }
    }
  }
}

/**
 *   Updates the machine state. Returns true if it was updated.
 *
 */
bool MaterialReader::next_state() 
{
  SCOPELOG(1);

  smatch match;

  CIMap<State> nextState = {
    { "engine",                  State::INITIAL  },
    { "porothermoelastic",       State::POROTHERMOELASTIC     },
    { "creep",                   State::CREEP     }
  };

  // Resets the state
  if (regex_match(line, RE_EMPTY)) { 
    current_state = State::INITIAL;
    return true; 
  }

  // New section. Changes the state
  if ( ! regex_search(line, match, RE_SEC) ) return false;

  string sec = match[1];

  if ( ! nextState.count(sec) ) flog << "Cannot find next state for section key '" << sec << "' (line " << ln << ").";

  // Save
  current_state  = nextState[ sec ];

  if ( regex_search(line, match, RE_SEC_NAME) ) 
  {
    string sname = match[2];
    /**/
    if ( iequals( sec, "engine" ) ) config.engine = sname;

    /**/
    if ( current_state == State::CREEP ) 
    {
      if ( iequals( sname , "carter" ) ) current_state = State::CREEP_CARTER;
      else flog << "Unknown creep model '" << sname << "'";
    }
  }

  dlog(1) << "NEXT STATE: " << current_state;

  return true;

}

/**
 *
 */
void MaterialReader::creep_carter_state() 
{
  SCOPELOG(1);
  smatch match;
  string vname;
  if ( regex_search( line, match, RE_STR_STR_STR ) ) 
    reg_param_str( match[1], match[2], match[3], "creep_carter" );

  else if ( regex_search( line, match, RE_STR_NUM ) ) 
    reg_param_dbl( match[1], "CON", stod(match[2]), "creep_carter" );

  else if ( regex_search( line, match, RE_STR_STR_NUM ) ) 
    reg_param_dbl( match[1], match[2], stod(match[3]), "creep_carter" );

}

/**
 *
 */
void MaterialReader::porothermoelastic_state() 
{
  SCOPELOG(1);
  smatch match;
  string vname;
  if ( regex_search( line, match, RE_STR_STR_STR ) ) 
    reg_param_str( match[1], match[2], match[3], "porothermoelastic" );

  else if ( regex_search( line, match, RE_STR_NUM ) ) 
    reg_param_dbl( match[1], "CON", stod(match[2]), "porothermoelastic" );

  else if ( regex_search( line, match, RE_STR_STR_NUM ) ) 
    reg_param_dbl( match[1], match[2], stod(match[3]), "porothermoelastic" );

}


/**
 *
 */
void MaterialReader::reg_param_str( string vname, string type, string val, string context )
{
  if ( ! KNOWN_VAR_TYPES.count(type) ) flog << "Unkwown variable type at '" << filename << "' (" << ln << "): " << type;
  if ( iequals( type, "file" ) ) config.file_param( vname, context ) = val;
}

/**
 *
 */
void MaterialReader::reg_param_dbl( string vname, string type, double val, string context )
{
  if ( iequals( type, "con" ) )  config.con_param(vname, context) = val ;
}
