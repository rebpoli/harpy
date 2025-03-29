
#include "config/ModelConfig.h"
#include "config/reader/ModelReader.h"
#include <iomanip>

/**
 *  
 */
ModelConfig::ModelConfig( string model_dir_ ) :
          model_dir( model_dir_ ), model_file( model_dir + "/MODEL") 
{ ModelReader( *this ); }


/**
 *
 */
ostream& operator<<(ostream& os, const ModelConfig & m)
{
  os << endl;
  os << "ModelConfig: " << endl;
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

  return os;
}
