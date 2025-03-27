
#include "config/ModelConfig.h"
#include "util/File.h"
#include <iomanip>
#include <regex>
#include "util/String.h"


/**
 *
 *
 */
ModelConfig::ModelConfig( string sys_name_ ) :
          sys_name( sys_name_ ), dir( sys_name ),
          model_file( dir + "/MODEL") 
{
  check_files();
  read_model();
}


/**
 *
 */
void ModelConfig::check_files()
{
  if ( ! dir_exists( dir ) )
    flog << "Cannot find model directory '" << dir << "'!";

  if ( ! file_exists( model_file ) )
    flog << "Cannot find model file '" << model_file << "'!";

  dlog(1) << "Model directory and file ok! We are good to go! (" << dir << ", " << model_file << ")";
}

/**
 *
 */
void ModelConfig::read_model()
{
  enum class State { INITIAL, SUBDOMAIN };
  ifstream file(model_file);
  if ( !file.is_open() ) flog << "Cannot open file " << model_file << ". Should have been resolved in check_files. What's wrong?";

  State currentState = State::INITIAL;

  // Regular expressions
  string tok = R"(([a-zA-Z_]+))";
  regex emptyRE ( R"(^\s*$)"               );
  regex keyRE   ( R"(\.)" + tok            );
  regex valueRE ( tok + R"(\s+)" + tok     );

  std::string line;
  while (std::getline(file, line)) 
  {
    line = remove_comments_and_trim( line );
    dlog(1) << line;

    // Empty line resets the machine
    if (regex_match(line, emptyRE)) 
    {
      dlog(1) << "  empty.";
      currentState = State::INITIAL;
      continue;
    }

    // Force state change if a known key is seen
    smatch match;
    if (regex_search(line, match, keyRE)) {
      string key = match[1];
      dlog(1) << "   key:" << key;
      if ( key == "SUBDOMAIN" ) 
      {
        currentState = State::SUBDOMAIN;
        continue;
      }
    }

    // Machine
    switch (currentState) 
    {
      //
      case State::INITIAL: 
      {
      }
      //
      case State::SUBDOMAIN: 
      {
        smatch match;
        if (regex_search(line, match, valueRE)) 
        {
          subdomain_material[ match[1] ] = match[2];
          continue;
        }
      }
      //
      default: break;
    } // switch
     
  } // machine loop
}

/**
 *
 */
ostream& operator<<(ostream& os, const ModelConfig & m)
{
  os << endl;
  os << "ModelConfig: (" << m.sys_name << ") : " << endl;
  os << "      " << setw(20) << left << "dir"  << ": " << m.dir << endl;
  os << "      " << setw(20) << left << "model file" << ": " << m.model_file << endl;
  os << "      " << "Materials:" << endl;
  for ( auto [ s, mat ] : m.subdomain_material )
    os << "           " << setw(15) << left << s << ": " << mat << endl;

  return os;
}
