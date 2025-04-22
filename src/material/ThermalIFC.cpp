
#include "material/ThermalIFC.h"

#include <iomanip>

/**
 *
 */
void ThermalIFC::reinit( uint eid, uint nqp ) 
{
  by_qp = &( by_elem[eid] );

  if ( ! by_qp->size() ) 
    by_qp->resize(nqp); 
}

/**
 *
 */
ostream& operator<<(ostream& os, const ThermalIFC & m)
{
  os << "THERMAL INTERFACE:" << endl;

  for ( auto & [ eid, pvec ] : m.by_elem )
  {
    os << "    EID:" << eid;
    os << "          [" << endl;
    for ( auto & p : pvec ) {
      os << "         temperature:" << setw(15) << p.temperature << endl;
    }
    os << "          ]" << endl;
  }

  return os;
}
