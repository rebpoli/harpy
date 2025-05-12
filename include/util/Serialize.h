
#pragma once

#include "libmesh/tensor_value.h"
#include "libmesh/vector_value.h"
#include "libmesh/point.h"
#include <boost/serialization/serialization.hpp>

#include <boost/serialization/vector.hpp>
#include <boost/serialization/map.hpp>
#include <boost/mpi.hpp>
namespace mpi = boost::mpi;

/**
 *
 *  Add serialization capabilities to a number of classes
 *
 */

namespace boost { 
namespace serialization {

/** **/
template<class Archive, typename T>
void serialize(Archive & ar, libMesh::TensorValue<T> & tensor, const unsigned int /*version*/ )
{
  for (unsigned int i = 0; i < 3; ++i) 
  for (unsigned int j = 0; j < 3; ++j) 
    ar & tensor(i,j);
} 
/** **/
template<class Archive, typename T>
void serialize(Archive & ar, libMesh::VectorValue<T> & vector, const unsigned int /*version*/ )
{ for (uint i=0; i<3; i++ ) ar & vector(i); }

/** **/
template<class Archive>
void serialize(Archive & ar, libMesh::Point & p, const unsigned int /*version*/ )
{ for (uint i=0; i<3; i++ ) ar & p(i); }


} } // namespace boost::serialization
