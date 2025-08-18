
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
struct MpiFileOps
{
  boost::mpi::communicator comm;
  MapT & the_map;

  MpiFileOps( MapT & a_map ) : comm(), the_map(a_map) {}

  // Tags
  enum { TAG_MAP = 100, TAG_OK = 101, TAG_KEYS = 102 };

  // API
  void save(const string& path) const;
  void load(const string& path);

private:
  friend class ser::access;

  template<class Ar>
  void serialize(Ar& ar, const unsigned /*version*/)
  { ar & ser::base_object<MapT>(the_map); }

  void write_file(const string& path) const;
  void read_file(const string& path) const;
};

/**
 *
 */
template<class MapT>
void MpiFileOps<MapT>::save(const string& path) const
{
  const int rank   = comm.rank();
  const int nprocs = comm.size();

  // serial
  if (nprocs == 1) { write_file(path); return; }

  ///
  if (rank == 0)
  {
    for (int src = 1; src < nprocs; ++src)
    {
      MapT shard;
      comm.recv(src, TAG_MAP, shard);
      the_map.insert(shard.begin(), shard.end());
    }

    write_file(path);

    int ok = 0;
    for (int dst = 1; dst < nprocs; ++dst)
      comm.send(dst, TAG_OK, ok);

    if (!ok) flog << "MapStore::save: root reported send failure.";
  }

  else // Rank != 0
  {
    comm.send(0, TAG_MAP, the_map);

    int ok = 0;
    comm.recv(0, TAG_OK, ok);

    if (!ok) flog << "MapStore::save: root reported write failure.";
  }
}

/**
 *
 */
template<class MapT>
void MpiFileOps<MapT>::load(const string& path)
{
  const int rank   = comm.rank();
  const int nprocs = comm.size();

  if (nprocs == 1) { read_file(path); return; }

  ///
  if (rank == 0)
  {
    read_file(path);
    for (int dst = 1; dst < nprocs; ++dst)
      comm.send(dst, TAG_MAP, the_map);
  }
  else // Rank != 0
  {
    the_map.clear();
    comm.recv(0, TAG_MAP, the_map);
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
  oa << the_map;
}

/**
 *
 */
template<class MapT>
void MpiFileOps<MapT>::read_file(const string& path) const
{
  ifstream ifs(path, ios::binary);
  if (!ifs) flog << "Cannot open '" << path << "'.";

  arch::binary_iarchive ia(ifs);
  ia >> the_map;
}

