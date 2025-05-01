
#include "material/ViscoplasticIFC.h"
#include "util/OutputOperators.h"

ViscoplasticIFC::~ViscoplasticIFC()
{
  for ( auto & [ _, vec ] : probes_by_pname )
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
ViscoplasticIFC::PropsTranspose::PropsTranspose( vector<Props> * by_qp )
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
  Props props;
  props.init_from_config (config, pt );

  ProbeIFC * probe_ifc = new ProbeIFC() ;
  probe_ifc->elem_id = eid;
  probe_ifc->pt = pt;
  probe_ifc->props = props;

  vector<ProbeIFC *> & vec = probes_by_pname[name];
  vec.push_back( probe_ifc );

  auto & vec_e = probes_by_elem[eid];
  vec_e.push_back( probe_ifc );   ///PROBE
};

/**
 *
 *
 */
void ViscoplasticIFC::Props::init_from_config( const MaterialConfig & config, const Point & pt )
{
  alpha_d          = config.get_property( "alpha_d",         pt,     "porothermoelastic" );
  beta_e           = config.get_property( "beta_e",          pt,     "porothermoelastic" );
  lame_mu          = config.get_property( "lame_mu",         pt,     "porothermoelastic" );
  lame_lambda      = config.get_property( "lame_lambda",     pt,     "porothermoelastic" );
  creep_carter_a   = config.get_property( "a",               pt,     "creep_carter" );
  creep_carter_q   = config.get_property( "q",               pt,     "creep_carter" );
  creep_carter_n   = config.get_property( "n",               pt,     "creep_carter" );

  sigtot = RealTensor();
  plastic_strain_n = RealTensor();
  plastic_strain   = RealTensor();
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
ostream& operator<<(ostream& os, const ViscoplasticIFC::Props & p)
{
  os << "VISCOPLASTIC PROPERTIES:" << endl;

  os << "    Static variables:" << endl;
  os << right << setw(20) << "lame_mu:" << setw(15) << p.lame_mu << endl;
  os << right << setw(20) << "lame_lambda:" << setw(15) << p.lame_lambda << endl;
  os << right << setw(20) << "alpha_d:" << setw(15) << p.alpha_d << endl;
  os << right << setw(20) << "beta_e:" << setw(15) << p.beta_e << endl;
  os << right << setw(20) << "creep_carter_a:" << setw(15) << p.creep_carter_a << endl;
  os << right << setw(20) << "creep_carter_q:" << setw(15) << p.creep_carter_q << endl;
  os << right << setw(20) << "creep_carter_n:" << setw(15) << p.creep_carter_n << endl;
  os << "    State variables:" << endl;
  os << right << setw(20) << "U:" << setw(15) << Print(p.U) << endl;
  os << right << setw(20) << "Temperature:" << setw(15) << p.temperature << endl;

  return os;
}
