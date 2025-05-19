
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
    if ( iequals( vname, "biot" ) )          return biot_file;
    if ( iequals( vname, "young" ) )         return young_file;
    if ( iequals( vname, "poisson" ) )       return poisson_file;
    if ( iequals( vname, "beta_e" ) )        return beta_e_file;
    if ( iequals( vname, "beta_d" ) )        return beta_d_file;
  }

  // Various modes of Munson Dawson
  if ( context == "creep_md1" ) { /* tbd */ }
  if ( context == "creep_md2" ) { /* tbd */ }

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
  if ( context == "creep_md1" )
  {
    if ( iequals( vname, "q" ) ) return creep_md1_q;
    if ( iequals( vname, "n" ) ) return creep_md1_n;
    if ( iequals( vname, "eps0" ) ) return creep_md1_eps0;
    if ( iequals( vname, "sig0" ) ) return creep_md1_sig0;

    if ( iequals( vname, "c" ) ) return creep_md1_c;
    if ( iequals( vname, "k" ) ) return creep_md1_k;
    if ( iequals( vname, "m" ) ) return creep_md1_m;
    if ( iequals( vname, "alpha_w" ) ) return creep_md1_alpha_w;
    if ( iequals( vname, "alpha_r" ) ) return creep_md1_alpha_r;
    if ( iequals( vname, "beta_w" ) ) return creep_md1_beta_w;
    if ( iequals( vname, "beta_r" ) ) return creep_md1_beta_r;
  }
  if ( context == "creep_md2" )
  {
    if ( iequals( vname, "q" ) ) return creep_md2_q;
    if ( iequals( vname, "n" ) ) return creep_md2_n;
    if ( iequals( vname, "eps0" ) ) return creep_md2_eps0;
    if ( iequals( vname, "sig0" ) ) return creep_md2_sig0;
    if ( iequals( vname, "c" ) ) return creep_md2_c;
    if ( iequals( vname, "k" ) ) return creep_md2_k;
    if ( iequals( vname, "m" ) ) return creep_md2_m;
    if ( iequals( vname, "alpha_w" ) ) return creep_md2_alpha_w;
    if ( iequals( vname, "alpha_r" ) ) return creep_md2_alpha_r;
    if ( iequals( vname, "beta_w" ) ) return creep_md2_beta_w;
    if ( iequals( vname, "beta_r" ) ) return creep_md2_beta_r;
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
  os << endl;
  os << "                            MaterialConfig for '" << m.name  << "'/ '" << m.cfg << "'" << endl;
  os << "                                 Filename '" << m.filename << "'" << endl;

  os << "                                 Thermoporoelastic:" << endl ;
  os << "                                    Por:              " << setw(15) << m.porosity     << setw(15) << m.porosity_file << endl;
  os << "                                    Permeability:     " << setw(15) << m.permeability << setw(15) << m.permeability_file << endl;
  os << "                                    Biot:             " << setw(15) << m.biot         << setw(15) << m.biot_file << endl;
  os << "                                    Young Modulus:    " << setw(15) << m.young        << setw(15) << m.young_file << endl;
  os << "                                    Poisson Coef:     " << setw(15) << m.poisson      << setw(15) << m.poisson_file << endl;
  os << "                                    Lame_mu:          " << setw(15) << m.lame_mu << endl;
  os << "                                    Lame_lambda:      " << setw(15) << m.lame_lambda << endl;
  os << "                                 Creep Model (MD1):" << endl ;
  os << "                                    EPS0:             " << setw(15) << m.creep_md1_eps0 << endl;
  os << "                                    SIG0:             " << setw(15) << m.creep_md1_sig0 << endl;
  os << "                                    Q:                " << setw(15) << m.creep_md1_q << endl;
  os << "                                    N:                " << setw(15) << m.creep_md1_n << endl;
  os << "                                    C:                " << setw(15) << m.creep_md1_c << endl;
  os << "                                    K:                " << setw(15) << m.creep_md1_k << endl;
  os << "                                    M:                " << setw(15) << m.creep_md1_m << endl;
  os << "                                    ALPHA_W:          " << setw(15) << m.creep_md1_alpha_w << endl;
  os << "                                    BETA_W:           " << setw(15) << m.creep_md1_beta_w << endl;
  os << "                                    ALPHA_R:          " << setw(15) << m.creep_md1_alpha_r << endl;
  os << "                                    BETA_R:           " << setw(15) << m.creep_md1_beta_r << endl;
  os << "                                 Creep Model (MD2):" << endl ;
  os << "                                    EPS0:             " << setw(15) << m.creep_md2_eps0 << endl;
  os << "                                    SIG0:             " << setw(15) << m.creep_md2_sig0 << endl;
  os << "                                    Q:                " << setw(15) << m.creep_md2_q << endl;
  os << "                                    N:                " << setw(15) << m.creep_md2_n << endl;
  os << "                                    C:                " << setw(15) << m.creep_md2_c << endl;
  os << "                                    K:                " << setw(15) << m.creep_md2_k << endl;
  os << "                                    M:                " << setw(15) << m.creep_md2_m << endl;
  os << "                                    ALPHA_W:          " << setw(15) << m.creep_md2_alpha_w << endl;
  os << "                                    BETA_W:           " << setw(15) << m.creep_md2_beta_w << endl;
  os << "                                    ALPHA_R:          " << setw(15) << m.creep_md2_alpha_r << endl;
  os << "                                    BETA_R:           " << setw(15) << m.creep_md2_beta_r << endl;

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
