
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
