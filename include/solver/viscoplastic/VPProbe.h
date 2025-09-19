#pragma once

#include "harpy/Global.h"
#include "solver/viscoplastic/VPProps.h"
#include <util/Serialize.h>
#include <set>
#include <boost/mpi.hpp>

namespace solver {
namespace viscoplastic {

namespace mpi = boost::mpi;

using libMesh::Point;

/// Probe interface
struct ProbeIFC { 
  ProbeIFC() = default;
  ProbeIFC( const ProbeIFC & ) = delete;
  uint elem_id; Point pt; VPProps props; uint pt_idx;
};

// Add serialization capabilities and MPI support
using ProbeByElemBaseMap = map< uint, vector<ProbeIFC *> >;
struct ProbeByElemMap : public ProbeByElemBaseMap
{
  boost::mpi::communicator world;
  void localize_to_one(ProbeByElemMap & global_map) ;

private:
  friend class boost::serialization::access;

  void _send();

  // On rank 0
  void _receive( ProbeByElemMap & global_map ) ;

  // For boost
  template<class Archive>
  void serialize(Archive& ar, const unsigned int /*ver*/) {
    ar & boost::serialization::base_object<ProbeByElemBaseMap>(*this);
  }
};

//
struct ProbeByPnameByElemMap : public map<string, ProbeByElemMap>
{
  boost::mpi::communicator world;

  vector<ProbeIFC *> probes_by_elem( uint eid )
  {
    vector<ProbeIFC *> ret;
    for ( auto & [ pname, m1 ] : *this )
      if  ( auto it = m1.find(eid); it != m1.end() )
        for ( auto & p : it->second ) 
          ret.push_back( p );
    return ret;
  }

  void localize_to_one(ProbeByPnameByElemMap& global_map)
  {
    uint rank   = world.rank();
    uint nprocs = world.size();
    if (nprocs == 1) { global_map = *this; return; }

    // 1) Collect local keys and send to root
    vector<string> local_keys;
    local_keys.reserve(this->size());
    for (auto const& kv : *this) local_keys.push_back(kv.first);

    if (rank == 0) 
    {
      // gather all key lists
      vector<vector<string>> all_keys;
      gather(world, local_keys, all_keys, /*root=*/0);

      // union + sort
      set<string> union_set;
      for (auto& v : all_keys) union_set.insert(v.begin(), v.end());
      vector<string> all_unique(union_set.begin(), union_set.end());

      // broadcast the canonical ordered list to everyone
      broadcast(world, all_unique, /*root=*/0);

      // 2) For each key, delegate to inner localize_to_one
      global_map.clear();
      for (auto const& key : all_unique) {
        // local copy if present; otherwise an empty map
        ProbeByElemMap local_inner;
        auto it = this->find(key);
        if (it != this->end()) local_inner = it->second;

        ProbeByElemMap merged_inner;
        local_inner.localize_to_one(merged_inner); // this performs MPI gather for this key

        // only root receives the merged map; store it
        global_map.emplace(key, move(merged_inner));
      }
    } else {
      // non-root: send keys to root
      gather(world, local_keys, /*root=*/0);

      // receive the canonical ordered list
      vector<string> all_unique;
      broadcast(world, all_unique, /*root=*/0);

      // For each key in the broadcast order, call inner localize_to_one in the same sequence
      for (auto const& key : all_unique) {
        ProbeByElemMap local_inner;
        auto it = this->find(key);
        if (it != this->end()) local_inner = it->second;

        ProbeByElemMap dummy; // ignored on non-root
        local_inner.localize_to_one(dummy);
      }
      // non-root: global_map remains untouched
    }
  }
};

}} // ns

// Add serialization for the specific types
namespace boost {
namespace serialization {
  using solver::viscoplastic::ProbeIFC;

  /** **/
  template<class Archive>
  void serialize(Archive & ar, ProbeIFC & p, const unsigned int /*ver*/)
  { 
    ar & p.elem_id & p.pt & p.props & p.pt_idx;  
  } 
  /** **/
}}  // Namespaces
