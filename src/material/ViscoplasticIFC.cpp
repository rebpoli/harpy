
#include "material/ViscoplasticIFC.h"

/**
 *
 */
void ViscoplasticIFC::reinit( uint eid, uint nqp )
{
  by_qp = &( by_elem[eid] );
  if ( ! by_qp->size() ) {
    by_qp->resize(nqp); valid=0; 
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

  os << right << setw(20) << "lame_mu:" << setw(15) << p.lame_mu << endl;
  os << right << setw(20) << "lame_lambda:" << setw(15) << p.lame_lambda << endl;
  os << right << setw(20) << "alpha_d:" << setw(15) << p.alpha_d << endl;
  os << right << setw(20) << "beta_e:" << setw(15) << p.beta_e << endl;
  os << right << setw(20) << "creep_carter_a:" << setw(15) << p.creep_carter_a << endl;
  os << right << setw(20) << "creep_carter_q:" << setw(15) << p.creep_carter_q << endl;
  os << right << setw(20) << "creep_carter_n:" << setw(15) << p.creep_carter_n << endl;

  return os;
}
