
#include "harpy/Coupler.h"
#include "util/OutputOperators.h"

/**
 *
 */
ostream& operator<<(ostream& os, const Coupler & m)
{
  os << endl;
  os << "Coupler:";
  for ( auto & [ eid, ecoup ] : m )
  {
    os << "    Element ID: '" << eid << "': " << endl;
    os << ecoup;
  }
  return os;
}

/**
 *
 */
ostream& operator<<(ostream& os, const ElemCoupler & m)
{
  os << "         (ElemCoupler) " << endl;
  for ( auto & [ vname, parbyqp ] : m.dbl_params)
    os << "             Var '" << vname << "': " << parbyqp << endl;

  return os;
}
