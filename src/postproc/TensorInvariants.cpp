#include "postproc/TensorInvariants.h"

/**
 * This object works on a DenseMatrix.
 * This constructor converts the RealTensor in a DenseMatrix
 * before initialization.
 */
TensorInvariants::TensorInvariants( const RealTensor & T ) :
  invarP(0), invarQ(0)
{
  DenseMatrix<Real> M(3,3);
  for ( uint i=0; i<3 ; i ++ )
  for ( uint j=0; j<3 ; j ++ )
    M(i,j)=T(i,j);

  init( M );
}

/**
 *  
 */
TensorInvariants::TensorInvariants( const DenseMatrix<Real> & M_ ) :
  invarP(0), invarQ(0)
{ 
  init(M_); 
}

/**
 *   Do all calculations of the invariants. 
 *   All the processing is here, then is accessed through getters.
 */
void TensorInvariants::init( const DenseMatrix<Real> & M_ )
{
  // Compute the eigenvalues/vectors
  DenseVector<Real> eval, eval_im;
  DenseMatrix<Real> evec;
  {
    DenseMatrix<Real> M = M_;
    M.evd_right( eval, eval_im, evec ); // Note: this func may change M
  }

  // Create the EigenPair
  for ( uint i = 0 ; i<3 ; i++ )
  {
    EigenPair ep;
    for ( uint j = 0 ; j<3 ; j++ ) ep.evec(j) = evec(j,i);

    ep.eval = eval(i);
    eigenPairs.push_back(ep);
  }

  // Sort  - the most compressive is the S1 (more negative)
  sort( eigenPairs.begin(), eigenPairs.end() );

  // P = 1/3 tr(\sigma)
  for ( uint i = 0 ; i < 3 ; i++ ) 
    invarP += M_(i, i) / 3.;

  // \sigma' = \sigma - p I
  // J2 = 1/2 tr( \sigma'_ij^2 )
  Real J2=0;

  DenseMatrix<Real> Sij = M_;
  for ( uint i = 0 ; i < 3 ; i++ ) 
    Sij(i,i) -= invarP;

  for ( uint i = 0 ; i < 3 ; i++ ) 
  for ( uint j = 0 ; j < 3 ; j++ ) 
    J2 += Sij(i,j) * Sij(i,j) / 2;

  invarQ = sqrt( J2 );
}

/**
 *
 *
 *
 */
ostream& operator<<(ostream& os, vector<TensorInvariants::EigenPair> & m)
{
  for ( auto ep : m ) os << ep.eval << " => " << ep.evec << endl;
  return os;
}



