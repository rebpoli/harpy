
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
    deviatoric.push_back( p.deviatoric );
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
