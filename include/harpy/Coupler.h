#pragma once

#include "base/Global.h"
#include "util/String.h"
#include <vector>
#include <map>

using harpy_string::CIMap;

/**
 *
 * This class holds the datastructures to couple two solvers.
 *
 * The structures are filled by the _update_coupler_ functions in the solvers.
 *
 */

class Solver;


/**
 *   The part of the coupler relevant to an element.
 *   To be used in the Material class.
 *
 */
struct ElemCoupler 
{
  template < typename T > struct Params : vector< T > { };

  CIMap<Params<double> > dbl_params;
  CIMap<Params<double> > dbl_params_old;  // Data from previous timestep

  ElemCoupler( uint eid_ ) : eid(eid_) {} 
  uint eid;
};


/**
 *      Map of ElemCoupler ( elem->id  :   ElemCoupler ) 
 */
struct  Coupler : map< uint , ElemCoupler >
{
  Coupler() {}
};

ostream& operator<<(ostream& os, const Coupler & m);
ostream& operator<<(ostream& os, const ElemCoupler & m);
