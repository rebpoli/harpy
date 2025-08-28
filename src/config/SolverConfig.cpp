
#include "config/SolverConfig.h"
#include "config/reader/SolverReader.h"
#include "util/OutputOperators.h"

#include <iomanip>

namespace config
{

/**
 *  
 */
SolverConfig::SolverConfig( string model_dir_, string sys_name_, string sys_cfg_ ) :
          model_dir( model_dir_ ), sys_file( model_dir + "/" + sys_name_),
          sys_name(sys_name_), sys_cfg(sys_cfg_)
{ 
  using namespace config::reader;
  SolverReader( *this ); 
}

/**
 *
 *
 */
bool SolverConfig::ExternalFile::check()
{ 
  if ( is_defined() ) return true;
  else flog << "External file is not completely defined!" << endl << *this;
  return false;
}

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

  for ( auto & [ v, fem ] : m.fem_by_var ) 
  {
    os << "    FEM (var: " << v << "):" << endl ;
    os << "       type:              " << setw(15) << fem.type     << endl;
    os << "       family:            " << setw(15) << fem.family   << endl;
    os << "       order:             " << setw(15) << fem.order    << endl;
    os << "       implicit:          " << setw(15) << fem.implicit     << endl;
  }

  os << m.external_file;

  return os;
}

/** **/
ostream& operator<<(ostream& os, const SolverConfig::ExternalFile & m)
{
  using util::operator<<;
  os << "    External_file: " << endl;
  os << "            filename:      " << setw(15) << m.filename << endl;
  os << "            grid_type:     " << setw(15) << m.grid_type << endl;
  os << "            grid_origin:   " << setw(15) << m.grid_origin << endl;
  os << "            min_radius:    " << setw(15) << m.min_radius << endl;
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

} // ns
