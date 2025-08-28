
#pragma once

#include "harpy/Global.h"

#include <map>
#include <vector>

namespace solver {
namespace thermal {

/** 
 *
 *    Holds the data to feed temperature into the material 
 *
 **/
struct ThermalIFC
{
  struct Props 
  { 
    double temperature; 
  };

  /// Data storage. By element, By Qp
  map< uint, vector<Props> > by_elem;

  /// Dynamic, set after initialization
  vector<Props> * by_qp;

  // Helpers
  void reinit( uint eid, uint nqp );

  // Getters
  Props & get( uint qp ) { return (*by_qp)[qp]; }
};


ostream& operator<<(ostream& os, const ThermalIFC & m);

}} // ns
