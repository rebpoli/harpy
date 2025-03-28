
#include "config/MaterialConfig.h"
#include "config/ReaderRegex.h"
#include "util/OutputOperators.h"

#include "util/File.h"
#include "util/String.h"
#include <iomanip>
#include <regex>
#include <set>


using namespace harpy_string;
using namespace MRDEF;

set<string> KNOWN_VAR_TYPES = { "CON", "FILE" };

/**
 *  
 */
MaterialConfig::MaterialConfig( const string & model_dir_, const string & name_, const string & cfg_ )   :
            model_dir(model_dir_), name(name_), cfg(cfg_) 
{
  filename = model_dir + "/" + name + "/" + to_upper_copy(cfg);

  check_files();
  parse_material_file();
}


/**
 *
 */
void MaterialConfig::check_files()
{
  if ( ! dir_exists( model_dir ) )
    flog << "Cannot find model directory '" << model_dir << "'!";

  if ( ! file_exists( filename ) )
    flog << "Cannot find material file '" << filename << "'!";

  dlog(1) << "Model directory and material file ok! We are good to go! (" << model_dir << ", " << filename << ")";
}

/**
 *  Reads the MATERIAL file. No need for a state machine as this is a linear file.
 */
void MaterialConfig::parse_material_file()
{
  ifstream file(filename);
  if ( !file.is_open() ) flog << "Cannot open file " << filename << ". Should have been resolved in check_files. What's wrong?";

  ln = 0;

  smatch match;
  while (  getline(file, line)  ) 
  {
    ln++;
    line = remove_comments_and_trim( line );
    if (regex_match(line, emptyRE)) continue;

    string vname;
    if ( regex_search( line, match, RE_STR_STR_STR ) ) 
      reg_param_str( match[1], match[2], match[3] );

    else if ( regex_search( line, match, RE_STR_NUM ) ) 
      reg_param_dbl( match[1], "CON", stod(match[2]) );

    else if ( regex_search( line, match, RE_STR_STR_NUM ) ) 
      reg_param_dbl( match[1], match[2], stod(match[3]) );
     
  } // Line loop
}


/**
 *
 */
void MaterialConfig::reg_param_str( string vname, string type, string val )
{
  if ( ! KNOWN_VAR_TYPES.count(type) ) flog << "Unkwown variable type at '" << filename << "' (" << ln << "): " << type;
  if ( iequals( type, "file" ) ) _file_param( vname ) = val;
}

/**
 *
 */
void MaterialConfig::reg_param_dbl( string vname, string type, double val )
{
  if ( iequals( type, "con" ) )  _con_param(vname) = val ;
}

/**
 *
 */
optional<string> & MaterialConfig::_file_param( string & vname )
{
  dlog(1) << "File param: " << vname;
  if ( iequals( vname, "porosity" ) )      return porosity_file;
  if ( iequals( vname, "permeability" ) )  return permeability_file;
  if ( iequals( vname, "biot" ) )          return biot_file;
  if ( iequals( vname, "bulk" ) )          return bulk_file;
  if ( iequals( vname, "skempton" ) )      return skempton_file;

  flog << "Variable name in Material '" << name << "', " << filename << " (" << ln << ")";
  return porosity_file;
}
/**
 *
 */
optional<double> & MaterialConfig::_con_param( string & vname )
{
  if ( iequals( vname, "porosity" ) )      return porosity;
  if ( iequals( vname, "permeability" ) )  return permeability;
  if ( iequals( vname, "biot" ) )          return biot;
  if ( iequals( vname, "bulk" ) )          return bulk;
  if ( iequals( vname, "skempton" ) )      return skempton;

  flog << "Variable name in Material '" << name << "', " << filename << " (" << ln << ")";
  return porosity;
}

/**
 *
 *
 */
ostream& operator<<(ostream& os, const MaterialConfig & m)
{
  os << endl;
  os << "                            MaterialConfig for '" << m.name  << "'/ '" << m.cfg << "'" << endl;
  os << "                                 Filename '" << m.filename << "'" << endl;

  os << "                                 Poroelastic:" << endl ;
  os << "                                    Por:              " << setw(15) << m.porosity     << setw(15) << m.porosity_file << endl;
  os << "                                    Permeability:     " << setw(15) << m.permeability << setw(15) << m.permeability_file << endl;
  os << "                                    Biot:             " << setw(15) << m.biot         << setw(15) << m.biot_file << endl;
  os << "                                    Bulk:             " << setw(15) << m.bulk         << setw(15) << m.bulk_file << endl;
  os << "                                    Skempton:         " << setw(15) << m.skempton     << setw(15) << m.skempton_file << endl;

  return os;
}
