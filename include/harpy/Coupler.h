#pragma once

#include "base/Global.h"
#include "util/String.h"
#include <vector>
#include <map>

#include "libmesh/vector_value.h"
#include "libmesh/tensor_value.h"

using namespace libMesh;

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
 *   A bunch of vectors indexed by the qp of the FE realization.
 *
 */
struct ElemCoupler 
{
  CIMap<vector<RealVectorValue>> vector_params;
  CIMap<vector<RealTensor> > tensor_params;

  CIMap<vector<double> > dbl_params;
  CIMap<vector<double> > dbl_params_old;  // Data from previous timestep

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
