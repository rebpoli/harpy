#pragma once

#include "harpy/Global.h"

#include "libmesh/transient_system.h"
#include "solver/viscoplastic/EGDataStruct.h"

namespace solver {
namespace viscoplastic {

class ViscoplasticSolver;
using namespace libMesh;

/**
 *
 */
class VPMatEG 
{
public:
  VPMatEG( TransientNonlinearImplicitSystem & sys_,
           ViscoplasticSolver & vpsolver_ );

  void init();                             /// Create datastructures from the mesh
                                            
  void residual_and_jacobian ( const NumericVector<Number> & soln, 
                               NumericVector<Number> * residual,
                               SparseMatrix<Number> * jacobian );

private:
  vector< EGFacePair > gamma_I;        /// Internal skeleton
  vector< EGFace > gamma_H;            /// External boundaries

  TransientNonlinearImplicitSystem & system;
  ViscoplasticSolver & vpsolver; 

  EGFEM fem_p, fem_n;                      /// Shape functions to be reinit'ed
  QGauss qrule;                            /// QRULE in the fem_p struct

  friend ostream& operator<<(ostream& os, const VPMatEG & m);
};

ostream& operator<<(ostream& os, const VPMatEG & m);

}} // ns
