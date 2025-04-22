
#include "material/ThermalIFC.h"


/**
 *
 */
void ThermalIFC::reinit( uint eid, uint nqp ) 
{
  by_qp = &( by_elem[eid] );

  if ( ! by_qp->size() ) 
    by_qp->resize(nqp); 
}
