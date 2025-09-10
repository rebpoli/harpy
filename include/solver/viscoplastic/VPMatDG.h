#pragma once

#include "harpy/Global.h"

#include "libmesh/transient_system.h"
#include "solver/viscoplastic/DGDataStruct.h"
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
class VPMatDG 
{
public:
  VPMatDG( TransientNonlinearImplicitSystem & sys_,
           ViscoplasticSolver & vpsolver_ );

  void init();                             /// Create datastructures from the mesh
                                            
  void residual_and_jacobian ( const NumericVector<Number> & soln, 
                               NumericVector<Number> * residual,
                               SparseMatrix<Number> * jacobian );

  void init_properties();
  void init_properties_gammad();

  // For Gamma_I (internal faces, skeleton)
  void reinit_I( DGFacePair & fp );
  void reinit_I( const NumericVector<Number> & soln , DGFacePair & fp );
  AD::Vec residual_gamma_I_qp( const AD::Vec & );
  void residual_and_jacobian_gammaI_qp( DGFacePair & fp );

  // For Gamma_D (external faces, Dirichlet)
  void reinit_D( DGFace & dgf );
  void reinit_D( const NumericVector<Number> & soln , DGFace & dgf );
  AD::Vec residual_gammaD_qp( const AD::Vec & );
  void residual_and_jacobian_gammaD_qp( DGFace & dgf );

  void update_gammad();

private:
  void setup_dofs( DGFace & dgf );
  void setup_dofs( DGFacePair & fp );
  inline bool next_qp();

  vector< DGFacePair > gamma_I;         /// Internal skeleton
  vector< DGDirichlet > gamma_D;

  TransientNonlinearImplicitSystem & system;
  ViscoplasticSolver & vpsolver; 

  DGFEM fem_p, fem_n;                      /// Shape functions to be reinit'ed
  QGauss qrule;                            /// QRULE in the fem_p struct
  VPProps * Pp, * Pn;                        /// Pointer to the resolved properties (see next_qp)
  vector<optional<double>> u_hat;

  /// This must include both elements (P and N concatenated)
  vector<dof_id_type> dof_indices;
  uint n_uvars, n_dofs_u;

  inline uint Ndof(uint i) { return n_dofs_u; }

  DenseMatrix<Number> Ke; /// Jacobian for the element
  DenseVector<Number> Re; /// RHS vector for the element

  uint QP;

  util::AD::ContextDG ad;

  // The size penalty
  double elem_penalty, beta;

  /** e: element ; i:dimension (x,y,z) ; B:element DOF**/
  inline uint idx_cg( uint e, uint i, uint B ) 
  { return e*3*n_dofs_cg + i*n_dofs_cg + B; }

  friend ostream& operator<<(ostream& os, const VPMatDG & m);
};

ostream& operator<<(ostream& os, const VPMatDG & m);

/** **/
inline bool VPMatDG::next_qp()
{ 
  if ( ++QP == qrule.n_points() ) return false;
  return true; 
}

}} // ns
