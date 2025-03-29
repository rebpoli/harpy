
#include "config/SystemConfig.h"
#include "config/reader/SystemReader.h"
#include <iomanip>

/**
 *  
 */
SystemConfig::SystemConfig( string model_dir_, string sys_name_, string sys_cfg_ ) :
          model_dir( model_dir_ ), sys_file( model_dir + "/" + sys_name_),
          sys_name(sys_name_), sys_cfg(sys_cfg_)
{ SystemReader( *this ); }


/**
 *
 */
ostream& operator<<(ostream& os, const SystemConfig & m)
{
  os << endl;
  os << "SystemConfig: " << endl;
  os << "    Material CFG_ID: " << endl;

  for ( auto [ s, mconf ] : m.material_config )
    os << "           " << setw(15) << left << s << ": " << mconf << endl;

  return os;
}

/**   **/
ostream& operator<<(ostream& os, const SystemConfig::MatConfigMap & m)
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
ostream& operator<<(ostream& os, const SystemConfig::MatConfig & m)
{
  os << setw(20) << m.name;
  os << setw(20) << m.cfg;
  return os;
}
