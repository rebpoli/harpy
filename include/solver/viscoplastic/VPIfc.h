#pragma once

#include "harpy/Global.h"
#include "config/MaterialConfig.h"
#include "solver/viscoplastic/VPProbe.h"
#include "solver/viscoplastic/VPProps.h"
#include "util/MpiFileOps.h"
#include <boost/serialization/access.hpp>

namespace solver {
namespace viscoplastic {

using namespace libMesh;

struct ViscoplasticIFC;

/// Full data
struct VPPropsByElemMap : public MpiFileOps< map<uint, vector<VPProps> > >
{ };

/** 
 * Holds the data to feed viscoplastic properties into/from the material 
 **/
struct ViscoplasticIFC
{
  ~ViscoplasticIFC();

  /// Data storage. By element, By Qp
  VPPropsByElemMap by_elem; 

  /// Dynamic, set after initialization
  vector<VPProps> * by_qp;

  /// 1 if the interface has been initialized
  bool valid;

  void reinit( uint eid, uint nqp=0 );

  /// Getters
  VPProps & get( uint qp ) { return (*by_qp)[qp]; }

  /// Probes organized by elemetns
  ProbeByPnameByElemMap probes_by_pname_by_elem;  // This map is MPI aware

  void add_probe_point( const MaterialConfig & config, string & name, uint eid, const Point & pt, uint pt_idx );

private:

  /* Serialization routines - polymorphic serialization. */
  friend class boost::serialization::access;
  template<class Ar> void serialize(Ar& ar, const unsigned /*version*/) 
  { }
};

/** OUTPUT STREAMS **/
ostream& operator<<(ostream& os, const ViscoplasticIFC & m);

}} // ns
