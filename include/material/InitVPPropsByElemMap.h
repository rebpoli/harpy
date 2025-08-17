#pragma once

#include "base/Global.h"
#include "util/Serialize.h"
#include "util/MpiFileOps.h"
#include "libmesh/tensor_value.h"

/** Snapshot of data - only the part that is needed from 
 * gravity activation (a stress initialization run) to the
 * analysis run.
 */
struct InitVPProps 
{ libMesh::RealTensor initial_strain; };

/** Serialization rule for InitVPProps **/
namespace boost { namespace serialization {

  template<class Archive>
  void serialize(Archive & ar, InitVPProps & p, const unsigned int /*ver*/)
  { ar & p.initial_strain; } 

} }  // Namespaces

/**
 * Serializeable map
 */

using InitVPPropsByElemMap_basemap = map< uint, vector<InitVPProps> >;
class InitVPPropsByElemMap : public InitVPPropsByElemMap_basemap
{ 
  MpiFileOps<InitVPPropsByElemMap_basemap> eng;

public:
  InitVPPropsByElemMap() : eng(*this) {}

  // Delegate to MpiFileOps
  void save( string fn ) { eng.save(fn); }
  void load( string fn ) { eng.load(fn); }
};

