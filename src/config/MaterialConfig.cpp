
#include "config/MaterialConfig.h"

#include "config/reader/MaterialReader.h"
#include "util/OutputOperators.h"

#include "util/String.h"


using harpy_string::iequals;
using harpy_string::to_upper_copy;
using harpy_string::ci_cmp;

/**
 *  
 */
MaterialConfig::MaterialConfig( const string & model_dir_, const string & name_, const string & cfg_ )   :
            model_dir(model_dir_), name(name_), cfg(cfg_) , filename( model_dir + "/" + name + "/" + to_upper_copy(cfg) )
{
  MaterialReader( *this );
}

/**
 *
 */
optional<string> & MaterialConfig::file_param( string & vname )
{
  dlog(1) << "File param: " << vname;
  if ( iequals( vname, "porosity" ) )      return porosity_file;
  if ( iequals( vname, "permeability" ) )  return permeability_file;
  if ( iequals( vname, "biot" ) )          return biot_file;
  if ( iequals( vname, "bulk" ) )          return bulk_file;
  if ( iequals( vname, "skempton" ) )      return skempton_file;

  flog << "Variable name in Material '" << name << "', " << filename;
  return porosity_file;
}
/**
 *
 */
optional<double> & MaterialConfig::con_param( string & vname )
{
  if ( iequals( vname, "porosity" ) )      return porosity;
  if ( iequals( vname, "permeability" ) )  return permeability;
  if ( iequals( vname, "biot" ) )          return biot;
  if ( iequals( vname, "bulk" ) )          return bulk;
  if ( iequals( vname, "skempton" ) )      return skempton;

  flog << "Variable name in Material '" << name << "', " << filename;
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

  for ( auto & [ v, fem ] : m.fem_by_var ) 
  {
    os << "                                 FEM (var: " << v << "):" << endl ;
    os << "                                    type:              " << setw(15) << fem.type     << endl;
    os << "                                    family:            " << setw(15) << fem.family   << endl;
    os << "                                    order:             " << setw(15) << fem.order    << endl;
  }

  return os;
}


/**
 * Enable the object to index a map.
 * The MaterialConfig is indexed by name and configuration
 */
bool MaterialConfig::operator<(const MaterialConfig & other) const
{
  ci_cmp c;
  if ( c( name, other.name ) ) return true;
  if ( c( other.name, name ) ) return false;
  return c( cfg, other.cfg );
}
