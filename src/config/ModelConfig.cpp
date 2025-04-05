
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
SolverConfig * ModelConfig::solver_config( string name )
{
  if ( ! solvers.count(name) ) 
  {
    string str; for ( auto & [ k, v ] : solvers ) { str += k + ", "; }
    flog << "Cannot find solver named '" << name << "'. Known solvers: " << str;
  }
  SolverConfig & sc = solvers.at(name);
  return &sc;
}

/**
 *
 */
ostream& operator<<(ostream& os, const CIMap<MaterialConfig> & m)
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
  for ( auto [ s, c ] : m.solvers )
    os << "           " << setw(15) << left << s << ": " << c << endl;

  os << "      " << "MATERIAL / Config:" << endl;
  for ( auto c : m.materials )
    os << "           " << setw(15) << left << c << endl;

  os << "      " << "Boundary  Config:" << endl;
  os << m.boundary_config;

  return os;
}
