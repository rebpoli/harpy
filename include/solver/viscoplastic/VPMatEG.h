#pragma once

#include "harpy/Global.h"

#include "libmesh/transient_system.h"
#include "solver/viscoplastic/EGDataStruct.h"
#include "util/Autodiff.h"

namespace solver {
namespace viscoplastic {

class ViscoplasticSolver;
using namespace libMesh;

//using namespace util;

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
  void residual_and_jacobian_qp();

  void reinit( const NumericVector<Number> & soln , EGFacePair & fp );

private:
  void setup_dofs( EGFacePair & fp );
  inline bool next_qp();

  vector< EGFacePair > gamma_I;        /// Internal skeleton
  vector< EGFace > gamma_H;            /// External boundaries

  TransientNonlinearImplicitSystem & system;
  ViscoplasticSolver & vpsolver; 

  EGFEM fem_p, fem_n;                      /// Shape functions to be reinit'ed
  QGauss qrule;                            /// QRULE in the fem_p struct

  /// This must include both elements (P and N concatenated)
  vector<dof_id_type> dof_indices_eg;
  uint n_dofsv, n_uvars, n_dofs_eg;

  /* For the internal loops:
   *       0,1,2 = Ux,Uy,Uz ; 3,4,5:UegX,UegY,UegZ
   */
  inline uint Ndof(uint i) { return i>2 ?  n_dofs_eg : n_dofsv; }

  DenseMatrix<Number> Ke; /// Jacobian for the element
  DenseVector<Number> Re; /// RHS vector for the element

  uint QP;
  util::AD::ContextEG ad;

  friend ostream& operator<<(ostream& os, const VPMatEG & m);
};

ostream& operator<<(ostream& os, const VPMatEG & m);

/** **/
inline bool VPMatEG::next_qp()
{ 
  if ( ++QP == qrule.n_points() ) return false;
  return true; 
}

}} // ns
