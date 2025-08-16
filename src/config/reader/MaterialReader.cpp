
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
      case State::CREEP_MD: { creep_md_state(); break; }
      case State::INITIALIZE: { initialize_state(); break; }
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
  smatch match;

  CIMap<State> nextState = {
    { "engine",                  State::INITIAL  },
    { "porothermoelastic",       State::POROTHERMOELASTIC     },
    { "creep",                   State::CREEP     },
    { "initialize",              State::INITIALIZE     }
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
      if      ( iequals( sname , "md" ) ) current_state = State::CREEP_MD;
      else                                 flog << "Unknown creep model '" << sname << "'";
    }
  }

  dlog(1) << "NEXT STATE: " << current_state;

  return true;

}

/**
 *
 */
void MaterialReader::initialize_state() 
{
  SCOPELOG(1);
  smatch match;

  if ( ! regex_search( line, match, RE_STR_STR ) ) return;
  string var_type = match[1]; // stress, pressure etc
  string method = match[2];

//  if (! config.creep_md ) config.creep_md = CreepMD{};
//  auto & creep_md = config.creep_md;

  // Steady state mechanism
  if ( iequals( var_type, "STRESS" ) ) 
  {
    using enum MatInitializeMethod;
    if ( iequals( method, "hydrostatic" ) ) config.initialize.method = HYDROSTATIC;
  }

}

/**
 *
 */
void MaterialReader::creep_md_state() 
{
  SCOPELOG(1);
  smatch match;

  // Must have something like SS [#] [#] [#] ...
  if ( ! regex_search( line, match, RE_STR_ETC ) ) return;
  string vname = match[1];

  if (! config.creep_md ) config.creep_md = CreepMD{};
  auto & creep_md = config.creep_md;

  // Steady state mechanism
  if ( iequals( vname, "SS" ) ) 
  {
    if ( ! regex_search( line, match, RE_CREEP_MD_SS ) ) flog << "Creep steady state line bad format (line " << ln << ")";
    double sig0     = stod( match[2] );
    double n        = stod( match[3] );
    double q        = stod( match[4] );
    double stretch  = 1;
    if ( match[5].length() ) stretch = stod( match[5] );
    creep_md->ss.push_back({sig0, n, q, stretch});
  }

  // Steady state mechanism
  if ( iequals( vname, "TR" ) ) 
  {
    if ( ! regex_search( line, match, RE_CREEP_MD_SS ) ) flog << "Creep steady state line bad format (line " << ln << ")";
    double sig0     = stod( match[2] );
    double m        = stod( match[3] );
    double c        = stod( match[4] );
    double alpha_w  = stod( match[5] );
    creep_md->tr.push_back({sig0, m, c, alpha_w});
  }
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
