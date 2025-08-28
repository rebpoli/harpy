
#include "config/MaterialConfig.h"
#include "config/reader/MaterialReader.h"
#include "util/OutputOperators.h"
#include "util/String.h"

namespace config
{

using util::iequals;
using util::to_upper_copy;
using util::ci_cmp;

/**
 *  
 */
MaterialConfig::MaterialConfig( const string & model_dir_, const string & name_, const string & cfg_ )   :
            model_dir(model_dir_), name(name_), cfg(cfg_) , filename( model_dir + "/" + name + "/" + to_upper_copy(cfg) )
{
  using namespace config::reader;
  MaterialReader( *this );
  
  // Compute the secondary properties
  if ( poisson && young )
  {
    lame_mu =  *young / 2 / ( 1 + *poisson );
    lame_lambda = *young * *poisson / (1 + *poisson) / (1 - 2 * (*poisson) );
    bulk_modulus = *young / 3 / ( 1 - 2* *poisson );
  }

  if ( bulk_modulus && beta_d ) alpha_d = *bulk_modulus * (*beta_d);
}

/**
 *
 */
optional<string> & MaterialConfig::file_param( string & vname, string context )
{
  dlog(1) << "File param: " << vname;
  if ( context == "porothermoelastic") 
  {
    if ( iequals( vname, "porosity" ) )      return porosity_file;
    if ( iequals( vname, "permeability" ) )  return permeability_file;
    if ( iequals( vname, "density" ) )       return density_file;
    if ( iequals( vname, "biot" ) )          return biot_file;
    if ( iequals( vname, "young" ) )         return young_file;
    if ( iequals( vname, "poisson" ) )       return poisson_file;
    if ( iequals( vname, "beta_e" ) )        return beta_e_file;
    if ( iequals( vname, "beta_d" ) )        return beta_d_file;
  }

  flog << "Variable name in Material '" << name << "', " << filename;
  return porosity_file;
}
/**
 *
 */
optional<double> & MaterialConfig::con_param( string & vname, string context ) 
{
  /** ** ** **/
  if ( context == "porothermoelastic" )
  {
    if ( iequals( vname, "porosity" ) )      return porosity;
    if ( iequals( vname, "permeability" ) )  return permeability;
    if ( iequals( vname, "biot" ) )          return biot;
    if ( iequals( vname, "density" ) )       return density;
    if ( iequals( vname, "young" ) )         return young;
    if ( iequals( vname, "poisson" ) )       return poisson;
    if ( iequals( vname, "beta_e" ) )        return beta_e;
    if ( iequals( vname, "beta_d" ) )        return beta_d;
    // secondary vars
    if ( iequals( vname, "lame_mu" ) )       return lame_mu;
    if ( iequals( vname, "lame_lambda" ) )   return lame_lambda;
    if ( iequals( vname, "bulk_modulus" ) )  return bulk_modulus;
    if ( iequals( vname, "alpha_d" ) )       return alpha_d;
  }

  /** ** ** **/
  flog << "Unknown variable name in Material '" << name << "': " << context << " . " << vname;
  return porosity;
} 
/** **/
const optional<double> & MaterialConfig::con_param( string & vname, string context ) const
{ return const_cast<optional<double>&>(const_cast<MaterialConfig*>(this)->con_param(vname,context)); }
const optional<string> & MaterialConfig::file_param( string & vname, string context ) const
{ return const_cast<optional<string>&>(const_cast<MaterialConfig*>(this)->file_param(vname,context)); }

/**
 *  TODO: We do not support file parameters yet. Only constants.
 */
void MaterialConfig::get_property( vector<double> & ret, string pname, const vector<Point> & xyz, string context ) const
{
  ret.clear();
  for ( auto & p : xyz ) 
    ret.push_back( get_property( pname, p, context ) );
}

/**
 *   TODO: support file and layer properties
 */
double MaterialConfig::get_property( string pname, const Point & xyz, string context ) const
{
  UNUSED(xyz);
  const optional<double> & prop = con_param(pname, context);
  if ( ! prop ) flog << "Property '" <<context <<"." << pname << "' is not defined for material '" << name << "'.";

  return *prop;
}

/**
 *
 *
 */
ostream& operator<<(ostream& os, const MaterialConfig & m)
{
  using util::operator<<;
  os << endl;
  os << "                            MaterialConfig for '" << m.name  << "'/ '" << m.cfg << "'" << endl;
  os << "                                 Filename '" << m.filename << "'" << endl;

  os << "                                 Thermoporoelastic:" << endl ;
  os << "                                    Por:              " << setw(15) << m.porosity     << setw(15) << m.porosity_file << endl;
  os << "                                    Permeability:     " << setw(15) << m.permeability << setw(15) << m.permeability_file << endl;
  os << "                                    Biot:             " << setw(15) << m.biot         << setw(15) << m.biot_file << endl;
  os << "                                    Density:          " << setw(15) << m.density      << setw(15) << m.density_file << endl;
  os << "                                    Young Modulus:    " << setw(15) << m.young        << setw(15) << m.young_file << endl;
  os << "                                    Poisson Coef:     " << setw(15) << m.poisson      << setw(15) << m.poisson_file << endl;
  os << "                                    Lame_mu:          " << setw(15) << m.lame_mu << endl;
  os << "                                    Lame_lambda:      " << setw(15) << m.lame_lambda << endl;
  os << "                                 Creep Model (MD):" << endl ;
  os << m.creep_md;
  return os;
}

ostream& operator<<(ostream& os, const CreepMD & m)
{
  for ( uint i=0 ; i < m.ss.size() ; i++ ) 
  {
    os << "                                   SteadyState (" << i << ")" << endl ;
    auto & ss = m.ss.at(i);
    os << "                                        sig0:     " << setw(15) << ss.sig0 << endl ;
    os << "                                        n:        " << setw(15) << ss.n << endl ;
    os << "                                        q:        " << setw(15) << ss.q << endl ;
    os << "                                        stretch:  " << setw(15) << ss.stretch << endl ;
  }

  for ( uint i=0 ; i < m.tr.size() ; i++ ) 
  {
    os << "                                   Transient (" << i << ")" << endl ;
    auto & tr = m.tr.at(i);
    os << "                                        sig0:     " << setw(15) << tr.sig0 << endl ;
    os << "                                        m:        " << setw(15) << tr.m << endl ;
    os << "                                        c:        " << setw(15) << tr.c << endl ;
    os << "                                        alpha_w:  " << setw(15) << tr.alpha_w << endl ;
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

} // ns
