#pragma once

#include "harpy/Global.h"

#include "libmesh/transient_system.h"
#include "solver/viscoplastic/EGDataStruct.h"
#include "util/Autodiff.h"

#include <optional>

namespace solver {
namespace viscoplastic {

class ViscoplasticSolver;
using namespace libMesh;

using namespace util;

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

  void init_properties();
  void init_properties_gammad();

  // For Gamma_I
  void reinit( EGFacePair & fp );
  void reinit( const NumericVector<Number> & soln , EGFacePair & fp );
  AD::Vec residual_qp( const AD::Vec & );
  void residual_and_jacobian_qp( EGFacePair & fp );

  // For Gamma_D
  void reinit( EGFace & egf );
  void reinit( const NumericVector<Number> & soln , EGFace & egf );
  AD::Vec residual_dirichlet_qp( const AD::Vec & );
  void residual_and_jacobian_dirichlet_qp( EGFace & egf );

  void update_gammad();

private:
  void setup_dofs( EGFace & egf );
  void setup_dofs( EGFacePair & fp );
  inline bool next_qp();

  vector< EGFacePair > gamma_I;         /// Internal skeleton
  vector< EGDirichlet > gamma_D;

  TransientNonlinearImplicitSystem & system;
  ViscoplasticSolver & vpsolver; 

  EGFEM fem_p, fem_n;                      /// Shape functions to be reinit'ed
  QGauss qrule;                            /// QRULE in the fem_p struct
  VPProps * Pp, * Pn;                        /// Pointer to the resolved properties (see next_qp)
  vector<optional<double>> u_hat;

  /// This must include both elements (P and N concatenated)
  vector<dof_id_type> dof_indices;
  uint n_uvars, n_dofs_cg, n_dofs_eg;

  /* For the internal loops:
   *       0,1,2 = Ux,Uy,Uz ; 3,4,5:UegX,UegY,UegZ
   */
  inline uint Ndof(uint i) { return i>2 ?  n_dofs_eg : n_dofs_cg; }

  DenseMatrix<Number> Ke; /// Jacobian for the element
  DenseVector<Number> Re; /// RHS vector for the element

  uint QP;

  // Autodiff holds only the EG part of the solution
  util::AD::ContextEG ad;

  // The size penalty
  double elem_penalty;

  /** e: element ; i:dimension (x,y,z) ; B:element DOF**/
  inline uint idx_cg( uint e, uint i, uint B ) 
  { return e*3*n_dofs_cg + i*n_dofs_cg + B; }

  friend ostream& operator<<(ostream& os, const VPMatEG & m);
};

ostream& operator<<(ostream& os, const VPMatEG & m);
ostream& operator<<(ostream& os, const EGDirichlet & m);

/** **/
inline bool VPMatEG::next_qp()
{ 
  if ( ++QP == qrule.n_points() ) return false;
  return true; 
}

}} // ns
