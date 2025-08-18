
#pragma once

#include "base/Global.h"

#include <map>
#include <set>
#include <vector>
#include <string>
#include <fstream>

#include <boost/serialization/access.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/set.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/mpi.hpp>


namespace mpi  = boost::mpi;
namespace ser  = boost::serialization;
namespace arch = boost::archive;

/**
 *
 * MpiFileOps: add save and load capabilities to a map of objects
 *             that is distributed in MPI processors.
 */

template<class MapT>
struct MpiFileOps : public MapT
{
  boost::mpi::communicator world;
  MpiFileOps() : world() {}

  // Tags
  enum { TAG_MAP = 100, TAG_OK = 101, TAG_KEYS = 102 };

  // API
  void save(const string& path);
  void load(const string& path);

  void localize_to_one( MapT & global_map ) const;
  void localize( MapT & global_map ) const;

private:
  friend class ser::access;

  template<class Ar>
  void serialize(Ar& ar, const unsigned /*version*/)
  { 
    MapT & self = static_cast<MapT&>(*this);
    ar & self; 
  }

  void write_file(const string& path) const;
  void read_file(const string& path);
};

/**
 *
 */
template<class MapT>
void MpiFileOps<MapT>::localize( MapT & global_map ) const
{
  MapT const& self = static_cast<MapT const&>(*this);

  const int rank   = world.rank();
  const int nprocs = world.size();

  // serial
  if (nprocs == 1) { global_map = self; return; }

  global_map.clear();

  std::vector<MapT> maps;

  all_gather(world, self, maps);

  global_map.clear();

  for (const auto& m : maps)
    for (const auto& [k, v] : m)
      global_map[k] = v;
}

/**
 *
 */
template<class MapT>
void MpiFileOps<MapT>::localize_to_one( MapT & global_map ) const
{
  MapT const& self = static_cast<MapT const&>(*this);

  const int rank   = world.rank();
  const int nprocs = world.size();

  // serial
  if (nprocs == 1) { global_map = self; return; }

  global_map.clear();

  // Sender
  if (rank != 0)
  { world.send(0,0, self); }

  // Receiver (root)
  else 
  {
    global_map = self; // Add my own data

    for ( uint i=1 ; i<nprocs ; ++i )
    {
      MapT remote_map;
      world.recv(i, 0, remote_map );

      for (const auto& [ key, val ] : remote_map) 
        global_map[key] = val;
    }
  }
}

/**
 *
 */
template<class MapT>
void MpiFileOps<MapT>::save(const string& path) 
{
  const int rank   = world.rank();
  const int nprocs = world.size();

  // serial
  if (nprocs == 1) { write_file(path); return; }

  ///
  if (rank == 0)
  {
    for (int src = 1; src < nprocs; ++src)
    {
      MapT shard;
      world.recv(src, TAG_MAP, shard);
      this->insert(shard.begin(), shard.end());
    }

    write_file(path);

    int ok = 0;
    for (int dst = 1; dst < nprocs; ++dst)
      world.send(dst, TAG_OK, ok);

    if (!ok) flog << "MapStore::save: root reported send failure.";
  }

  else // Rank != 0
  {
    world.send(0, TAG_MAP, *this);

    int ok = 0;
    world.recv(0, TAG_OK, ok);

    if (!ok) flog << "MapStore::save: root reported write failure.";
  }
}

/**
 *
 */
template<class MapT>
void MpiFileOps<MapT>::load(const string& path)
{
  const int rank   = world.rank();
  const int nprocs = world.size();

  if (nprocs == 1) { read_file(path); return; }

  ///
  if (rank == 0)
  {
    read_file(path);
    for (int dst = 1; dst < nprocs; ++dst)
      world.send(dst, TAG_MAP, *this);
  }
  else // Rank != 0
  {
    this->clear();
    world.recv(0, TAG_MAP, *this);
  }
}

/**
 *
 */
template<class MapT>
void MpiFileOps<MapT>::write_file(const string& path) const
{
  ofstream ofs(path, ios::binary);
  if (!ofs) flog << "Cannot open '" << path << "'.";

  arch::binary_oarchive oa(ofs);
  oa << *this;
}

/**
 *
 */
template<class MapT>
void MpiFileOps<MapT>::read_file(const string& path)
{
  ifstream ifs(path, ios::binary);
  if (!ifs) flog << "Cannot open '" << path << "'.";

  arch::binary_iarchive ia(ifs);
  ia >> *this;
}

