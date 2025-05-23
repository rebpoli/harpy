
#include "material/ViscoplasticIFC.h"
#include "util/OutputOperators.h"

ViscoplasticIFC::~ViscoplasticIFC()
{
  for ( auto & [ _, m1 ] : probes_by_pname_by_elem )
  for ( auto & [ _, vec ] : m1 )
  for ( auto p : vec ) delete(p);
}

/**
 *
 */
void ViscoplasticIFC::reinit( uint eid, uint nqp )
{
  by_qp = &( by_elem[eid] );

  if ( nqp )
  if ( ! by_qp->size() ) 
  {
    by_qp->resize(nqp); 
    valid=0; 
  }
}

/**
 *   Creates a transposed datastructure from a by_qp vector.
 */
PropsTranspose::PropsTranspose( vector<VPProps> * by_qp )
{
  for ( auto & p : *by_qp ) {
    sigtot.push_back( p.sigtot );
    sigeff.push_back( p.sigeff );
    von_mises.push_back( p.von_mises );
    epskk.push_back( p.epskk );
    F.push_back( p.F );
    deviatoric.push_back( p.deviatoric );
    plastic_strain.push_back( p.plastic_strain );
    plastic_strain_rate.push_back( p.plastic_strain_rate );
  }
}

/**
 *
 */
void ViscoplasticIFC::add_probe_point( const MaterialConfig & config, string & name, uint eid, const Point & pt )
{
  VPProps props;
  props.init_from_config (config, pt );

  ProbeIFC * probe_ifc = new ProbeIFC() ;
  probe_ifc->elem_id = eid;
  probe_ifc->pt = pt;
  probe_ifc->props = props;

  auto & m1 = probes_by_pname_by_elem[name];
  auto & vec_e = m1[eid];
  vec_e.push_back( probe_ifc );   ///PROBE
};

/**
 *
 *
 */
void VPProps::init_from_config( const MaterialConfig & config, const Point & pt )
{
  alpha_d            = config.get_property( "alpha_d",         pt,     "porothermoelastic" );
  beta_e             = config.get_property( "beta_e",          pt,     "porothermoelastic" );
  lame_mu            = config.get_property( "lame_mu",         pt,     "porothermoelastic" );
  lame_lambda        = config.get_property( "lame_lambda",     pt,     "porothermoelastic" );

  creep_md1_q        = config.get_property( "q",               pt,     "creep_md1" );
  creep_md1_n        = config.get_property( "n",               pt,     "creep_md1" );
  creep_md1_eps0     = config.get_property( "eps0",            pt,     "creep_md1" );
  creep_md1_sig0     = config.get_property( "sig0",            pt,     "creep_md1" );

  creep_md1_c        = config.get_property( "c",               pt,     "creep_md1" );
  creep_md1_k        = config.get_property( "k",               pt,     "creep_md1" );
  creep_md1_m        = config.get_property( "m",               pt,     "creep_md1" );
  creep_md1_alpha_w  = config.get_property( "alpha_w",         pt,     "creep_md1" );
  creep_md1_alpha_r  = config.get_property( "alpha_r",         pt,     "creep_md1" );
  creep_md1_beta_w   = config.get_property( "beta_w",          pt,     "creep_md1" );
  creep_md1_beta_r   = config.get_property( "beta_r",          pt,     "creep_md1" );

//  creep_md2_q      = config.get_property( "q",               pt,     "creep_md2" );
//  creep_md2_n      = config.get_property( "n",               pt,     "creep_md2" );
//  creep_md2_eps0   = config.get_property( "eps0",            pt,     "creep_md2" );
//  creep_md2_sig0   = config.get_property( "sig0",            pt,     "creep_md2" );

  sigtot = RealTensor();
  plastic_strain_n = RealTensor();
  plastic_strain   = RealTensor();
  creep_md1_etr   = 0;
  creep_md1_etr_n = 0;
}

/**
 * Updates the properties based on precalculated U and GRAD_U
 *
 */
void VPProps::update( const RealVectorValue & U_, const RealTensor & GRAD_U_,
                      double dt )
{
  U=U_; GRAD_U = GRAD_U_;

  RealTensor sigeff;
  for (uint i=0; i<3; i++) for (uint j=0; j<3; j++) 
  for (uint k=0; k<3; k++) for (uint l=0; l<3; l++)
    sigeff(i,j) += C_ijkl(i,j,k,l) * ( GRAD_U(k,l) - plastic_strain(k,l) );

  sigtot = sigeff;
  for (uint k=0; k<3; k++ ) 
    sigtot(k,k) -= alpha_d * ( temperature - initial_temperature );

  deviatoric = sigtot;
  for (uint i=0; i<3; i++ )
  for (uint k=0; k<3; k++ ) 
    deviatoric(i,i) -= (1./3.) * sigtot(k,k);

  double J2=0;
  for (uint i=0; i<3; i++ )
  for (uint j=0; j<3; j++ ) 
    J2 += (1./2.) * deviatoric(i,j) * deviatoric(i,j);

  von_mises = sqrt( 3 * J2 );

  // Compute plastic strain rate
  double R_ = 8.3144;   // Universal gas constant [ J/mol/K ]

  double etr_star = creep_md1_k * exp( creep_md1_c * temperature ) * 
                      pow( ( von_mises/creep_md1_sig0 ), creep_md1_m );
  double zeta = 0;
  if ( abs(etr_star) > 1e-20 ) 
    zeta = 1 - creep_md1_etr / etr_star;

  double alpha = creep_md1_alpha_w, beta = creep_md1_beta_w;
  if ( zeta <= 0 ) { alpha = creep_md1_alpha_r; beta = creep_md1_beta_r; }

  F = exp( alpha * zeta * zeta ) * pow( ( von_mises / creep_md1_sig0 ), beta * zeta * zeta );

  double etr_rate = ( F - 1 ) * creep_md1_eps0 *
                      exp( - creep_md1_q / R_ / temperature ) *
                      pow( von_mises/creep_md1_sig0 , creep_md1_n );
  creep_md1_etr = creep_md1_etr_n + etr_rate * dt;

  /// NOTE: This must match the calculations in the residual_qp function (ViscoPlasticMaterial)
  plastic_strain_rate =  3./2. * F * creep_md1_eps0 *
                         exp( - creep_md1_q / R_ / temperature ) *
                         pow( von_mises/creep_md1_sig0 , creep_md1_n-1 ) *
                         deviatoric * (1./creep_md1_sig0) ;

  plastic_strain = dt * plastic_strain_rate;
  for (uint i=0; i<3; i++ )
  for (uint j=0; j<3; j++ ) 
    plastic_strain(i,j) += plastic_strain_n(i,j);
}

/**
 *
 */
ostream& operator<<(ostream& os, const ViscoplasticIFC & m)
{
  os << "VISCOPLASTIC INTERFACE:" << endl;

  for ( auto & [ eid, pvec ] : m.by_elem )
  {
    os << "    EID:" << eid;
    os << "          [" << endl;
    for ( auto & p : pvec ) {
      os << "             lame_mu:" << setw(15) << p.lame_mu << endl;
      os << "         lame_lambda:" << setw(15) << p.lame_lambda << endl;
    }
    os << "          ]" << endl;
  }

  return os;
}

/**
 *
 */
ostream& operator<<(ostream& os, const VPProps & p)
{
  os << "VISCOPLASTIC PROPERTIES:" << endl;

  os << "    Static variables:" << endl;
  os << right << setw(20) << "lame_mu:" << setw(15) << p.lame_mu << endl;
  os << right << setw(20) << "lame_lambda:" << setw(15) << p.lame_lambda << endl;
  os << right << setw(20) << "alpha_d:" << setw(15) << p.alpha_d << endl;
  os << right << setw(20) << "beta_e:" << setw(15) << p.beta_e << endl;

  os << right << setw(20) << "creep_md1_eps0:" << setw(15) << p.creep_md1_eps0 << endl;
  os << right << setw(20) << "creep_md1_sig0:" << setw(15) << p.creep_md1_sig0 << endl;
  os << right << setw(20) << "creep_md1_q:" << setw(15) << p.creep_md1_q << endl;
  os << right << setw(20) << "creep_md1_n:" << setw(15) << p.creep_md1_n << endl;
  os << right << setw(20) << "creep_md1_c:" << setw(15) << p.creep_md1_c << endl;
  os << right << setw(20) << "creep_md1_k:" << setw(15) << p.creep_md1_k << endl;
  os << right << setw(20) << "creep_md1_m:" << setw(15) << p.creep_md1_m << endl;
  os << right << setw(20) << "creep_md1_alpha_w:" << setw(15) << p.creep_md1_alpha_w << endl;
  os << right << setw(20) << "creep_md1_beta_w:" << setw(15) << p.creep_md1_beta_w << endl;
  os << right << setw(20) << "creep_md1_alpha_r:" << setw(15) << p.creep_md1_alpha_r << endl;
  os << right << setw(20) << "creep_md1_beta_r:" << setw(15) << p.creep_md1_beta_r << endl;

  os << right << setw(20) << "creep_md2_eps0:" << setw(15) << p.creep_md2_eps0 << endl;
  os << right << setw(20) << "creep_md2_sig0:" << setw(15) << p.creep_md2_sig0 << endl;
  os << right << setw(20) << "creep_md2_q:" << setw(15) << p.creep_md2_q << endl;
  os << right << setw(20) << "creep_md2_n:" << setw(15) << p.creep_md2_n << endl;
  os << right << setw(20) << "creep_md2_c:" << setw(15) << p.creep_md2_c << endl;
  os << right << setw(20) << "creep_md2_k:" << setw(15) << p.creep_md2_k << endl;
  os << right << setw(20) << "creep_md2_m:" << setw(15) << p.creep_md2_m << endl;
  os << right << setw(20) << "creep_md2_alpha_w:" << setw(15) << p.creep_md2_alpha_w << endl;
  os << right << setw(20) << "creep_md2_beta_w:" << setw(15) << p.creep_md2_beta_w << endl;
  os << right << setw(20) << "creep_md2_alpha_r:" << setw(15) << p.creep_md2_alpha_r << endl;
  os << right << setw(20) << "creep_md2_beta_r:" << setw(15) << p.creep_md2_beta_r << endl;
  os << "    State variables:" << endl;
  os << right << setw(20) << "U:" << setw(15) << Print(p.U) << endl;
  os << right << setw(20) << "Temperature:" << setw(15) << p.temperature << endl;
  os << right << setw(20) << "Plastic_strain:" << setw(15) << Print(p.plastic_strain) << endl;
  os << right << setw(20) << "Plastic_strain_n:" << setw(15) << Print(p.plastic_strain_n) << endl;
  os << right << setw(20) << "Plastic_strain_rate:" << setw(15) << Print(p.plastic_strain_rate) << endl;
  return os;
}
