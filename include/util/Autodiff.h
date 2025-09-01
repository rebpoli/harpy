
#pragma once

#include "harpy/Global.h"

/*
 *     Sets a couple of aliases to make the code more readable.
 */

#include <autodiff/forward/real.hpp>
#include <autodiff/forward/real/eigen.hpp>
#include <autodiff/common/eigen.hpp>

#include "libmesh/tensor_value.h"

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
  struct ContextEG
  {

    /** n_elem_: 2 if it considers 2 neighbors **/
    ContextEG( uint n_uvars_ ) : n_dofs(0), n_dofsv(0), n_uvars(n_uvars_),  n_dofs_eg(0), n_elem(0)  {}

    /** **/
    inline void init(  uint n_dofsv_, uint n_dofs_eg_, uint n_elem_=1 ) 
    { 
      n_dofsv = n_dofsv_;
      n_dofs_eg = n_dofs_eg_;
      n_elem = n_elem_;

      // Number of DoFs per element
      n_dofs = n_dofsv * 3 + n_dofs_eg*(n_uvars-3);

      uint tot_dofs = n_dofs * n_elem;

      ad_Uib.resize( tot_dofs ); 
      ad_Uib.setZero();
      ad_Fib.resize( tot_dofs );
      ad_Fib.setZero();
      ad_Jijbm.resize( tot_dofs, tot_dofs );
      ad_Jijbm.setZero();
    }

    /** e: element ; i:dimension (x,y,z) ; B:element DOF**/
    inline uint idx( uint e, uint i, uint B ) 
    { 
      uint ret = e * n_dofs + B;

      if ( i > 2 ) 
        ret += 3*n_dofsv + (i-3) * n_dofs_eg;
      else 
        ret += i*n_dofsv ;

      return ret; 
    }

    /* Single element shortcuts */
    inline AD::real & Uib( uint i, uint B ) { return Ueib( 0, i, B ); }
    inline AD::real & Fib( uint i, uint B ) { return Feib( 0, i, B ); }
    inline AD::real Jijbm( uint i, uint j, uint B, uint M ) { return Jeijbm( 0, i, j, B, M ) ; }

    /* Multiple elements */
    inline AD::real & Ueib( uint e, uint i, uint B )
    { return ad_Uib[ idx(e, i,B) ]; }
    inline AD::real & Feib( uint e, uint i, uint B )
    { return ad_Fib[ idx(e, i,B) ]; }
    inline AD::real Jeijbm( uint e, uint i, uint j, uint B, uint M )
    { return ad_Jijbm( idx(e, i,B), idx(e, j,M) ); }

  private:

    /* dof counters */
    uint n_dofs, n_dofsv, n_uvars, n_dofs_eg, n_elem;

    /* Flatened stuff */
    AD::Vec ad_Uib, ad_Fib; 
    AD::Mat ad_Jijbm;
  };

} } // ns
