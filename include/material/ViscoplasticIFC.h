#pragma once

#include "base/Global.h"
#include "config/MaterialConfig.h"

#include "material/ProbeIFC.h"
#include "material/VPProps.h"
#include "material/InitVPPropsByElemMap.h"

#include <boost/serialization/access.hpp>

using namespace libMesh;
struct ViscoplasticIFC;

/// Full data
using VPPropsByElemMap = map< uint, vector<VPProps> >;

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

  void add_probe_point( const MaterialConfig & config, string & name, uint eid, const Point & pt );

  void save_initial_strain( const string & filename );
  void load_initial_strain( const string & filename );

private:
  InitVPPropsByElemMap snapshot_initial_strain();

  /* Serialization routines - polymorphic serialization. */
  friend class boost::serialization::access;
  template<class Ar> void serialize(Ar& ar, const unsigned /*version*/) 
  { }
};

/** OUTPUT STREAMS **/
ostream& operator<<(ostream& os, const ViscoplasticIFC & m);
