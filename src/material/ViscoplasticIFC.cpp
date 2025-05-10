
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
  alpha_d          = config.get_property( "alpha_d",         pt,     "porothermoelastic" );
  beta_e           = config.get_property( "beta_e",          pt,     "porothermoelastic" );
  lame_mu          = config.get_property( "lame_mu",         pt,     "porothermoelastic" );
  lame_lambda      = config.get_property( "lame_lambda",     pt,     "porothermoelastic" );

  creep_md1_q      = config.get_property( "q",               pt,     "creep_md1" );
  creep_md1_n      = config.get_property( "n",               pt,     "creep_md1" );
  creep_md1_eps0   = config.get_property( "eps0",            pt,     "creep_md1" );
  creep_md1_sig0   = config.get_property( "sig0",            pt,     "creep_md1" );

//  creep_md2_q      = config.get_property( "q",               pt,     "creep_md2" );
//  creep_md2_n      = config.get_property( "n",               pt,     "creep_md2" );
//  creep_md2_eps0   = config.get_property( "eps0",            pt,     "creep_md2" );
//  creep_md2_sig0   = config.get_property( "sig0",            pt,     "creep_md2" );

  sigtot = RealTensor();
  plastic_strain_n = RealTensor();
  plastic_strain   = RealTensor();
}

/**
 *
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
  //
  // Compute plastic strain rate
  double R_ = 8.3144;   // Universal gas constant [ J/mol/K ]

  plastic_strain_rate =  3./2. * creep_md1_eps0 *
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

  os << right << setw(20) << "creep_md2_eps0:" << setw(15) << p.creep_md2_eps0 << endl;
  os << right << setw(20) << "creep_md2_sig0:" << setw(15) << p.creep_md2_sig0 << endl;
  os << right << setw(20) << "creep_md2_q:" << setw(15) << p.creep_md2_q << endl;
  os << right << setw(20) << "creep_md2_n:" << setw(15) << p.creep_md2_n << endl;
  os << "    State variables:" << endl;
  os << right << setw(20) << "U:" << setw(15) << Print(p.U) << endl;
  os << right << setw(20) << "Temperature:" << setw(15) << p.temperature << endl;

  return os;
}
