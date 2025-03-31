
#include "config/ModelConfig.h"
#include "config/reader/ModelReader.h"
#include <iomanip>

/**
 *
 * Global instance of the model configuration, read in
 * load time.
 *
 */
ModelConfig * MODEL;

/**
 *  
 */
ModelConfig::ModelConfig( string model_dir_ ) :
          model_dir( model_dir_ ), model_file( model_dir + "/MODEL") 
{ ModelReader( *this ); }


/**
 *
 */
ostream& operator<<(ostream& os, const map<string,MaterialConfig> & m)
{
  os << "Map of MaterialConfig:" << endl;
  for ( auto & [ k , m ] : m )
    os << "     " << k << " => " << m;
  return os;
}

/**
 *
 */
ostream& operator<<(ostream& os, const ModelConfig & m)
{
  os << endl;
  os << "ModelConfig: " << endl;
  os << "      " << setw(20) << left << "systemloop"  << ": " << m.systemloop << endl;
  os << "      " << setw(20) << left << "model_dir"  << ": " << m.model_dir << endl;
  os << "      " << setw(20) << left << "model file" << ": " << m.model_file << endl;
  os << "      " << "Timestep:" << endl;
  for ( auto [ s, val ] : m.timestep )
    os << "           " << setw(15) << left << s << ": " << val << endl;
  os << "      " << "SYSTEM / Config:" << endl;
  for ( auto [ s, c ] : m.systems )
    os << "           " << setw(15) << left << s << ": " << c << endl;

  os << "      " << "MATERIAL / Config:" << endl;
  for ( auto [ s, c ] : m.materials )
    os << "           " << setw(15) << left << s << ": " << c << endl;

  os << "      " << "Boundary  Config:" << endl;
  os << m.boundary_config;

  return os;
}
