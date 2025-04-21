#pragma once

#include "base/Global.h"
#include "util/String.h"
#include <vector>
#include <map>

#include "libmesh/vector_value.h"
#include "libmesh/tensor_value.h"

using namespace libMesh;


///**  Templated generic property coupler with access functions **/
//template< typename T >
//struct PCoupler 
//{
//  PCoupler( uint nqp_ ) : nqp(nqp_ ) {}
//  map< uint, vector<T> > data;
//  vector<T> & get( uint eid ) { 
//    if ( ! data.count(eid) ) data[eid] = vector<T>(3);
//    return data.at[eid];
//  }
//  virtual string name() = 0;
//};

///**  Properties for coupling temperature (ex: Temperature from ThermalSolver to MechanicalSolver **/
//struct ThermalCouplerProperties 
//{
//  ThermalCouplerProperties() :      temperature(0) 
//  {}
//  double temperature; 
//};

///**
// * Specialization of the generic property coupler to thermal coupler
// * Prefering to use inheritance for better debugging. Could be a typedef.
// */
//struct ThermalCoupler : public PCoupler<ThermalCouplerProperties> 
//{ virtual string name() { return "ThermalCoupler"; } };



using harpy_string::CIMap;

/**
 *
 * This class holds the datastructures to couple two solvers.
 *
 * The structures are filled by the _update_coupler_ functions in the solvers.
 *
 */

class Solver;

///// DEPRECATING THESE //////


/**
 *   The part of the coupler relevant to an element.
 *   To be used in the Material class.
 *
 *   A bunch of vectors indexed by the qp of the FE realization.
 *
 */
//struct ElemCoupler 
//{
//  CIMap<vector<RealVectorValue>> vector_params;
//  CIMap<vector<RealTensor> > tensor_params;

//  CIMap<vector<double> > dbl_params;
//  CIMap<vector<double> > dbl_params_old;  // Data from previous timestep

//  ElemCoupler( uint eid_ ) : eid(eid_) {} 
//  uint eid;

//  void clear( vector<string> vnames ) {
//    for ( auto vname : vnames ) if ( dbl_params.count(vname) ) dbl_params.at(vname).clear();
//    for ( auto vname : vnames ) if ( tensor_params.count(vname) ) tensor_params.at(vname).clear();
//    for ( auto vname : vnames ) if ( vector_params.count(vname) ) vector_params.at(vname).clear();
//  }
//};


///**
// *      Map of ElemCoupler ( elem->id  :   ElemCoupler ) 
// */
//struct  Coupler : map< uint , ElemCoupler >
//{
//  Coupler() {}

//  ElemCoupler & elem_coupler( uint eid ) {
//    if ( ! count( eid ) ) flog << "Coupler has not been initialized for element '" << eid << "'. Something is terribly wrong.";
//    return at(eid);
//  }
//};

//ostream& operator<<(ostream& os, const vector<RealTensor> & m);
//ostream& operator<<(ostream& os, const Coupler & m);
//ostream& operator<<(ostream& os, const ElemCoupler & m);
