#pragma once

#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"
#include "libmesh/tensor_value.h"
#include "harpy/Global.h"
#include <vector>

namespace postproc {
namespace stress {

using namespace libMesh;

class TensorInvariants 
{
  private:
    // To help ordering the eigenPairs
    class EigenPair 
    {
      public:
        EigenPair() : evec(3), eval(0) {}
        DenseVector<Real> evec;
        Real eval;
        bool operator<(const EigenPair & other) { return eval < other.eval; }
    };

    // Ordered [0]:S1 ; [1]:S2 ; [2]: S3
    vector<EigenPair> eigenPairs;
    Real invarP, invarQ;

    // Initializes the internal structures (called from the constructor)
    void init( const DenseMatrix<Real> & M_ );
    
  public:

    TensorInvariants( const RealTensor & T ) ;
    TensorInvariants( const DenseMatrix<Real> & M_ ) ;

    // Convenient getters
    // i:0-2 => S1-S3    j:0-2 => X-Z
    double get_Sij( uint i, uint j ) { return eigenPairs[i].evec(j); }
    void get_Si( uint i, DenseVector<Real> & evec, Real & eval ) { evec = eigenPairs[i].evec; eval = eigenPairs[i].eval; }
    void get_S1( DenseVector<Real> & evec, Real & eval ) { get_Si( 0, evec, eval ); }
    void get_S2( DenseVector<Real> & evec, Real & eval ) { get_Si( 1, evec, eval ); }
    void get_S3( DenseVector<Real> & evec, Real & eval ) { get_Si( 2, evec, eval ); }

    double get_P() { return invarP; }
    double get_Q() { return invarQ; }

    // eval = eigen value
    double S1_eval() { return eigenPairs[0].eval; }
    double S2_eval() { return eigenPairs[1].eval; }
    double S3_eval() { return eigenPairs[2].eval; }

    // Convenience
    void get_eigenvectors( DenseVector<Real> & S1, DenseVector<Real> & S2, DenseVector<Real> & S3 ) { 
      S1 = eigenPairs[0].evec;  S2 = eigenPairs[1].evec;  S3 = eigenPairs[2].evec;  
    }

    friend ostream& operator<<(ostream& os, vector<EigenPair> & m);
};

ostream& operator<<(ostream& os, vector<TensorInvariants::EigenPair> & m);

}} // ns
