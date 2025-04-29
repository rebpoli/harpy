
#include "config/reader/InoutReader.h"

#include "config/reader/ReaderRegex.h"
#include "util/String.h"
#include "util/File.h"

using namespace MRDEF;
using namespace harpy_string;
using namespace libMesh;

/**
 *
 */
InoutReader::InoutReader( InoutConfig & _config ) :
  config(_config), filename(config.model_dir + "/INOUT"),
  curr_probe_config(0), ln(0), line(""), current_state(INITIAL)
{
  check_files();
  parse_inout_file();
}


/**
 *
 */
void InoutReader::check_files()
{
  string model_dir = config.model_dir;

  if ( ! dir_exists( model_dir ) ) flog << "Cannot find model directory '" << model_dir << "'!";

  if ( ! file_exists( filename ) ) flog << "Cannot find material file '" << filename << "'!";

  dlog(1) << "Model directory and material file ok! We are good to go! (" << model_dir << ", " << filename << ")";
}

/**
 *  Reads the INOUT file.
 */
void InoutReader::parse_inout_file()
{
  SCOPELOG(1);

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
      case State::PROBE: { probe_state(); break; }
      default: { wlog << "Ignoring line " << ln << " >> " << line; break; }
    }
  }
}

/**
 *
 */
bool InoutReader::next_state() 
{
  smatch match;

  CIMap<State> nextState = {
    { "probe",                  State::PROBE  }
  };

  // Resets the state
  if (regex_match(line, RE_EMPTY)) { 
    current_state = State::INITIAL;
    curr_probe_config = 0;
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
    if ( iequals( sec, "probe" ) ) 
      curr_probe_config = & ( config.probes.emplace_back(sname) );

  }

  dlog(1) << "NEXT STATE: " << current_state;

  return true;
}

/**
 *
 */
void InoutReader::probe_state() 
{
  SCOPELOG(1);
  smatch match;
  string vname;

  // Types
  CISet PT_PROPS    = { "from", "to", "center", "normal" };
  CISet UINT_PROPS  = { "npts" };
  CISet DBL_PROPS   = { "radius", "dtheta" };
  CISet STR_PROPS   = { "type", "boundary", "order" };

  if ( ! regex_search( line, match, RE_STR_ETC ) ) 
  { wlog << "Unmatched line? [" << line << "] (" << ln << ")"; return; }

  string k = match[1];

  /** Points **/
  if ( PT_PROPS.count(k) )
  {
    if ( ! regex_search( line, match, RE_STR_NUM_NUM_NUM) ) flog << "Wrong format at line " << ln << ": [" << line << "]";

    Point v = Point(stod(match[2]), stod(match[3]), stod(match[4]));
    if ( iequals( k, "from" )  )  curr_probe_config->from = v;
    if ( iequals( k, "to" ) )     curr_probe_config->to = v;
    if ( iequals( k, "center" ) ) curr_probe_config->center = v;
    if ( iequals( k, "normal" ) ) curr_probe_config->normal = v;

    return;
  }

  /** Unsigned ints **/
  if ( UINT_PROPS.count(k) )
  {
    if ( ! regex_search( line, match, RE_STR_UINT) ) flog << "Wrong format at line " << ln << ": [" << line << "]";

    uint v = stoi(match[2]);
    if ( iequals( k, "npts" ) ) curr_probe_config->npts = v;
    return;
  }

  /** Doubles **/
  if ( DBL_PROPS.count(k) )
  {
    if ( ! regex_search( line, match, RE_STR_NUM) ) flog << "Wrong format at line " << ln << ": [" << line << "]";
    dlog(1) << "match NUM - >" << line << "<";
     double v = stod(match[2]);
    if ( iequals( k, "radius" ) ) curr_probe_config->radius = v;
    if ( iequals( k, "dtheta" ) ) curr_probe_config->dtheta = v;
    return;
  }

  /** Doubles **/
  if ( STR_PROPS.count(k) )
  {
    if ( ! regex_search( line, match, RE_STR_STR) ) flog << "Wrong format at line " << ln << ": [" << line << "]";
    string v = match[2];
    if ( iequals( k, "type" ) )      curr_probe_config->type = v;
    if ( iequals( k, "boundary" ) )  curr_probe_config->boundaries.push_back(v);
    if ( iequals( k, "order" )  )    curr_probe_config->order = v;
    return;
  }
}
