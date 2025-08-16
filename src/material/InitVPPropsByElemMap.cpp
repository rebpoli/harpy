#include "material/InitVPPropsByElemMap.h"

#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/mpi/collectives.hpp>

#include <fstream>

/**
 *
 */
InitVPPropsByElemMap::InitVPPropsByElemMap() : 
  InitVPPropsByElemBaseMap(),
  world(),
  rank(world.rank()), 
  size(world.size()) 
{}

/**
 * This function is MPI aware
 */
void InitVPPropsByElemMap::write_to_file( const string & filename )
{
  InitVPPropsByElemMap global_map;
  localize_to_one( global_map );

  if ( rank == 0 )
  {
    ofstream ofs( filename, ios::binary );
    boost::archive::binary_oarchive oa( ofs );
    oa << global_map;
  }
}

/**
 * This function is MPI aware
 */
void InitVPPropsByElemMap::read_from_file( const string & filename )
{
  if ( rank == 0 ) 
  {
    ifstream ifs( filename, ios::binary );
    if ( !ifs ) flog << "Failed to open '" << filename << "'.";

    boost::archive::binary_iarchive ia( ifs );
    InitVPPropsByElemMap local_map;
    ia >> *this;
  }

  mpi::broadcast( world, *this, 0 );
}

/**
 * 
 **/
void InitVPPropsByElemMap::localize_to_one( InitVPPropsByElemMap & global_map )
{
  // single proc ?
  if (size == 1) 
  { 
    if (rank == 0) global_map = *this;
    return; 
  } 

  global_map.clear();
  if (rank != 0) 
    _send(); 
  else 
    _receive( global_map );

  world.barrier();   // sync
}

/**
 *
 */
void InitVPPropsByElemMap::_send() {
  world.send(0, 0, *this);  
  return; 
}

/**
 *
 */
void InitVPPropsByElemMap::_receive( InitVPPropsByElemMap & global_map ) 
{
  global_map = *this;  // add my own data

  for (int i = 1; i < size; ++i) 
  {
    InitVPPropsByElemMap remote_map;
    world.recv(i, 0, remote_map);

    for (const auto& [ key, vec ] : remote_map) 
    {
      auto & tvec = global_map[key];  // create if non existing
      tvec.insert(tvec.end(), vec.begin(), vec.end());
    }
  }
}


