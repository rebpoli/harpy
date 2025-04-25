
#pragma once

/*
 *     Sets a couple of aliases to make the code more readable.
 */

#include <autodiff/forward/real.hpp>
#include <autodiff/forward/real/eigen.hpp>
#include <autodiff/common/eigen.hpp>

#include "libmesh/tensor_value.h"

namespace AD
{
  using real = autodiff::real;
  using Mat = autodiff::MatrixXreal;
  using Vec = autodiff::VectorXreal;
  using autodiff::jacobian;

  inline void dump( const Mat & src, libMesh::RealTensor & trg )
  {
    for ( uint i=0; i<src.rows() ; i++ )
    for ( uint j=0; j<src.cols() ; j++ )
      trg(i,j) = val( src(i,j) );
  }
}
