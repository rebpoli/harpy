
#include "config/ModelConfig.h"
#include "util/File.h"
#include <iomanip>


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
}

/**
 *
 */
ostream& operator<<(ostream& os, const ModelConfig & m)
{
  os << "ModelConfig: (" << m.sys_name << ") : " << endl;
  os << "      " << setw(20) << left << "dir"  << ": " << m.dir << endl;
  os << "      " << setw(20) << left << "model file" << ": " << m.model_file << endl;
  return os;
}
