
#include "config/SolverConfig.h"
#include "config/reader/SolverReader.h"
#include <iomanip>

/**
 *  
 */
SolverConfig::SolverConfig( string model_dir_, string sys_name_, string sys_cfg_ ) :
          model_dir( model_dir_ ), sys_file( model_dir + "/" + sys_name_),
          sys_name(sys_name_), sys_cfg(sys_cfg_)
{ SolverReader( *this ); }


/**
 *
 */
ostream& operator<<(ostream& os, const SolverConfig & m)
{
  os << endl;
  os << "SolverConfig: " << endl;
  os << "    Model_dir       : " << m.model_dir << endl;
  os << "    Sys_file        : " << m.sys_file << endl;
  os << "    Sys_name        : " << m.sys_name << endl;
  os << "    Sys_cfg         : " << m.sys_cfg << endl;
  os << "    mesh_filename   : " << m.mesh_filename << endl;
  os << "    Material CFG_ID : " << endl;

  for ( auto [ s, mconf ] : m.mat_config_by_name )
    os << "           " << setw(15) << left << s << ": " << mconf << endl;

  return os;
}

/**   **/
ostream& operator<<(ostream& os, const SolverConfig::MatConfigMap & m)
{
  os << "------------------------------------------------------------------------------------------" << endl;
  os << setw(20) << "Subdomain";
  os << setw(20) << "Material";
  os << setw(20) << "Config";
  os << endl;
  for ( auto & [ sd, mconf ] : m )
  {
    os << setw(20) << sd;
    os <<  mconf << endl;
  }
  return os;
}
ostream& operator<<(ostream& os, const SolverConfig::MatNameAndCfg & m)
{
  os << setw(20) << m.name;
  os << setw(20) << m.cfg;
  return os;
}
