
#include "solver/viscoplastic/VPProbe.h"

namespace solver {
namespace viscoplastic {

/**
 *      SERIALIZATION STUFF
 */
void ProbeByElemMap::localize_to_one(ProbeByElemMap & global_map)
{
  uint rank = world.rank();
  if (world.size() == 1) { if (rank == 0) global_map = *this; return; }  // single proc.

  global_map.clear();
  if (rank != 0) _send(); else _receive( global_map );
  world.barrier();   // sync
}

/*
 * On rank != 0
 */
void ProbeByElemMap::_send() {
  world.send(0, 0, *this);  
  return; 
}

/*
 * On rank = 0
 */
void ProbeByElemMap::_receive( ProbeByElemMap & global_map ) 
{
  global_map = *this;  // add my own data

  for (int i = 1; i < world.size(); ++i) 
  {
    ProbeByElemMap remote_map;
    world.recv(i, 0, remote_map);

    for (const auto& [ key, vec ] : remote_map) 
    {
      auto & tvec = global_map[key];  // create if non existing
      tvec.insert(tvec.end(), vec.begin(), vec.end());
    }
  }
}

}} // ns
