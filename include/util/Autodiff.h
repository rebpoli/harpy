
#pragma once

#include "harpy/Global.h"

/*
 *     Sets a couple of aliases to make the code more readable.
 */

#include <autodiff/forward/real.hpp>
#include <autodiff/forward/real/eigen.hpp>
#include <autodiff/common/eigen.hpp>

#include "libmesh/tensor_value.h"
#include "libmesh/point.h"

namespace util {
namespace AD {
  /** Types **/
  using real = autodiff::real;
  using Mat = autodiff::MatrixXreal;
  using Vec = autodiff::VectorXreal;

  using autodiff::jacobian;

  /** **/
  inline void dump( const AD::Mat & src, libMesh::RealTensor & trg )
  {
    for ( uint i=0; i<src.rows() ; i++ )
    for ( uint j=0; j<src.cols() ; j++ )
      trg(i,j) = val( src(i,j) );
  }

  /** **/
  inline double norm( const AD::Mat & src )
  {
    libMesh::RealTensor trg;
    for ( uint i=0; i<src.rows() ; i++ )
    for ( uint j=0; j<src.cols() ; j++ )
      trg(i,j) = val( src(i,j) );
    return trg.norm();
  }

  /** **/
  inline AD::Vec dot( const AD::Mat & t, const AD::Vec & n )
  {
    AD::Vec ret = AD::Vec::Zero(n.size());
    for ( uint i=0; i<t.rows() ; i++ )
    for ( uint j=0; j<t.cols() ; j++ )
      ret(i) += t(i,j) * n(j);
    return ret;
  }

  /** **/
  inline AD::Vec dot( const AD::Mat & t, const libMesh::Point & n )
  {
    AD::Vec ret = AD::Vec::Zero(3);
    for ( uint i=0; i<t.rows() ; i++ )
    for ( uint j=0; j<t.cols() ; j++ )
      ret(i) += t(i,j) * n(j);
    return ret;
  }

  /** **/
  struct ContextEG
  {

    /** n_elem_: 2 if it considers 2 neighbors **/
    ContextEG( uint n_uvars_ ) : n_elem(0), n_dofs(0), n_uvars(n_uvars_), n_dofs_cg(0), n_dofs_eg(0)  {}

    /** **/
    inline void init(  uint n_dofs_cg_, uint n_dofs_eg_, uint n_elem_=1 ) 
    { 
      n_dofs_cg = n_dofs_cg_;
      n_dofs_eg = n_dofs_eg_;
      n_elem = n_elem_;

      // Number of DoFs per element
      n_dofs = n_dofs_cg * 3 + n_dofs_eg*(n_uvars-3);

      uint tot_dofs = n_dofs * n_elem;

      ad_Uib.resize( tot_dofs ); 
      ad_Uib.setZero();
      ad_Fib.resize( tot_dofs );
      ad_Fib.setZero();
      ad_Jijbm.resize( tot_dofs, tot_dofs );
      ad_Jijbm.setZero();
    }

    /** e: element ; i:dimension (x,y,z) ; B:element DOF**/
    inline uint idx( uint i, uint B ) { return idx(0,i,B); }
    inline uint idx( uint e, uint i, uint B ) 
    { 
      uint ret = e * n_dofs + B;

      if ( i > 2 ) 
        ret += 3*n_dofs_cg + (i-3) * n_dofs_eg;
      else 
        ret += i*n_dofs_cg ;

//      dlog(1) << "idx e(" << e << ") i(" << i << ") B(" << B << "): " << ret << " /// n_dofs(" << n_dofs << ") n_dofs_cg(" << n_dofs_cg << ") n_dofs_eg(" << n_dofs_eg << ")";
      return ret; 
    }

    /* Single element shortcuts */
    inline AD::real & Uib( uint i, uint B ) { return Ueib( 0, i, B ); }
    inline AD::real & Fib( uint i, uint B ) { return Feib( 0, i, B ); }
    inline AD::real Jijbm( uint i, uint j, uint B, uint M ) { return Jenijbm( 0, 0, i, j, B, M ) ; }

    /* Multiple elements */
    inline AD::real & Ueib( uint e, uint i, uint B )
    { return ad_Uib[ idx(e,i,B) ]; }
    inline AD::real & Feib( uint e, uint i, uint B )
    { return ad_Fib[ idx(e,i,B) ]; }
    // 
    inline AD::real Jenijbm( uint e, uint n, uint i, uint j, uint B, uint M )
    { return ad_Jijbm( idx(e,i,B), idx(n,j,M) ); }

    /* dof counters */
    uint n_elem, n_dofs, n_uvars, n_dofs_cg, n_dofs_eg;

    /* Flatened stuff */
    AD::Vec ad_Uib, ad_Fib; 
    AD::Mat ad_Jijbm;
  };

} } // ns
