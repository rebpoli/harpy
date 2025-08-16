#pragma once

#include "base/Global.h"
#include "material/VPProps.h"
#include <util/Serialize.h>

/// Probe interface
struct ProbeIFC { 
  ProbeIFC() = default;
  ProbeIFC( const ProbeIFC & ) = delete;
  uint elem_id; Point pt; VPProps props; 
};

// Add serialization capabilities and MPI support
using ProbeByElemBaseMap = map< uint, vector<ProbeIFC *> >;
struct ProbeByElemMap : public ProbeByElemBaseMap
{
  ProbeByElemMap() : ProbeByElemBaseMap(), world(), rank(world.rank()), size(world.size()) {}

  /// Main function that collects data from all processes and localizes to rank 0
  void localize_to_one(ProbeByElemMap & global_map) ;

private:
  mpi::communicator world;
  int rank, size;
  friend class boost::serialization::access;

  void _send();

  // On rank 0
  void _receive( ProbeByElemMap & global_map ) ;

  // For boost
  template<class Archive>
  void serialize(Archive& ar, const unsigned int version) {
    UNUSED(version);
    ar & boost::serialization::base_object<ProbeByElemBaseMap>(*this);
  }
};

using ProbeByPnameByElemMap = map< string, ProbeByElemMap > ;

// Add serialization for the specific types
namespace boost {
namespace serialization {
  /** **/
  template<class Archive>
  void serialize(Archive & ar, ProbeIFC & p, const unsigned int version)
  { 
    UNUSED(version);
    ar & p.elem_id; ar & p.pt; ar & p.props;  
  } 
  /** **/
} }  // Namespaces
