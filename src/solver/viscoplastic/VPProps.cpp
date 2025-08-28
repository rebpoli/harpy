#include "solver/viscoplastic/VPProps.h"
#include "util/OutputOperators.h"

namespace solver {
namespace viscoplastic {

using util::Print;

/**
 *
 *
 */
void VPProps::init_from_config( const MaterialConfig & config, const Point & pt )
{
  density            = config.get_property( "density",         pt,     "porothermoelastic" );
  biot               = config.get_property( "biot",            pt,     "porothermoelastic" );
  alpha_d            = config.get_property( "alpha_d",         pt,     "porothermoelastic" );
  beta_e             = config.get_property( "beta_e",          pt,     "porothermoelastic" );
  lame_mu            = config.get_property( "lame_mu",         pt,     "porothermoelastic" );
  lame_lambda        = config.get_property( "lame_lambda",     pt,     "porothermoelastic" );

  if ( config.creep_md ) creep_md = *(config.creep_md);

  sigtot           = RealTensor();
  plastic_strain_n = RealTensor();
  plastic_strain   = RealTensor();
  initial_stress   = RealTensor();
}

/**
 * Updates the properties based on precalculated U and GRAD_U
 *
 */
void VPProps::update( const RealVectorValue & U_, const RealTensor & GRAD_U_,
                      double dt )
{
  U=U_; GRAD_U = GRAD_U_;

  RealTensor sigeff = initial_stress;
  for (uint i=0; i<3; i++) for (uint j=0; j<3; j++) 
  for (uint k=0; k<3; k++) for (uint l=0; l<3; l++)
    sigeff(i,j) += 
      C_ijkl(i,j,k,l) *
      ( GRAD_U(k,l) - plastic_strain(k,l) );

  sigtot = sigeff;
  for (uint k=0; k<3; k++ ) 
  {
    sigtot(k,k) -= alpha_d * ( temperature - initial_temperature );
    sigtot(k,k) -= biot    * ( pressure - initial_pressure );
  }

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

  // STEADY STATE CREEP
  double creep_ss_rate = 0;
  for ( auto & ss : creep_md.ss ) 
    creep_ss_rate +=  exp( - pow(ss.q / R_ / temperature, ss.stretch) ) *
                      pow( von_mises/ss.sig0 , ss.n );
                 

  // TRANSIENT CREEP (find F)
  double F = 1;
  for ( auto & tr : creep_md.tr )
  {
    double etr_star = exp( tr.c * temperature ) * 
                      pow( ( von_mises/tr.sig0 ), tr.m );

    double zeta = 0;
    if ( abs(etr_star) > 1e-20 ) zeta = 1 - creep_md.etr / etr_star;

    double alpha = tr.alpha_w;
    if ( zeta <= 0 ) alpha = 0;

    F *= exp( alpha * zeta * zeta ) ;
  }

  double etr_rate = ( F - 1 ) * creep_ss_rate;
  creep_md.etr = creep_md.etr_n + etr_rate * dt;

  /// NOTE: This must match the calculations in the residual_qp function (ViscoPlasticMaterial)
  plastic_strain_rate = 0;
  if ( von_mises )
  for ( auto & ss : creep_md.ss ) 
    plastic_strain_rate += 
                      3./2. * F * 
                      exp( - pow(ss.q / R_ / temperature, ss.stretch) ) *
                      pow( von_mises/ss.sig0 , ss.n-1 ) * 
                      deviatoric / ss.sig0;

  plastic_strain = dt * plastic_strain_rate;
  for (uint i=0; i<3; i++ )
  for (uint j=0; j<3; j++ ) 
    plastic_strain(i,j) += plastic_strain_n(i,j);
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
    pressure.push_back( p.pressure );
    epskk.push_back( p.epskk );
    F.push_back( p.F );
    deviatoric.push_back( p.deviatoric );
    plastic_strain.push_back( p.plastic_strain );
    plastic_strain_rate.push_back( p.plastic_strain_rate );
    initial_stress.push_back( p.initial_stress );

    // Sign convention: compresion is negative
    RealTensor s_ = p.sigtot;
    for ( uint i=0; i<3; i++ ) s_(i,i) += p.pressure;
    sigeff_terz.push_back( s_ );
  }
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
  os << right << setw(20) << "biot:" << setw(15) << p.biot << endl;
  os << right << setw(20) << "alpha_d:" << setw(15) << p.alpha_d << endl;
  os << right << setw(20) << "beta_e:" << setw(15) << p.beta_e << endl;

  os << right << setw(20) << "creep_md:" << setw(15) << p.creep_md << endl;
  os << "    State variables:" << endl;
  os << right << setw(20) << "U:" << setw(15) << Print(p.U) << endl;
  os << right << setw(20) << "Initial Temperature:" << setw(15) << p.initial_temperature << endl;
  os << right << setw(20) << "Temperature:" << setw(15) << p.temperature << endl;
  os << right << setw(20) << "Initial Pressure:" << setw(15) << p.initial_pressure << endl;
  os << right << setw(20) << "Pressure:" << setw(15) << p.pressure << endl;
  os << right << setw(20) << "Plastic_strain:" << setw(15) << Print(p.plastic_strain) << endl;
  os << right << setw(20) << "Plastic_strain_n:" << setw(15) << Print(p.plastic_strain_n) << endl;
  os << right << setw(20) << "Plastic_strain_rate:" << setw(15) << Print(p.plastic_strain_rate) << endl;
  return os;
}

}} // ns
