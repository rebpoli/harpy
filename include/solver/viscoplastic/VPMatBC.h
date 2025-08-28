#pragma once

#include "harpy/Global.h"

#include "util/Autodiff.h"
#include <optional>
#include "libmesh/quadrature_gauss.h"
#include "libmesh/transient_system.h"

namespace config { class MaterialConfig; } 

namespace solver {
namespace viscoplastic {

class ViscoplasticSolver;
using config::MaterialConfig;
using namespace libMesh;

/**
 *
 *  THE MATERIAL BOUNDARY CLASS:
 *  Inherits as much as possible from the material.
 *  The FEM structures are one dimension lower
 *
 */
class ViscoPlasticMaterialBC 
{
public:
  ViscoPlasticMaterialBC( suint sid_, const MaterialConfig * config_,
                          TransientNonlinearImplicitSystem & sys_, 
                          ViscoplasticSolver & vpsolver_ );

  void reinit( const NumericVector<Number> & soln, const Elem & elem_, uint side );
  void residual_and_jacobian_qp();
  void residual_and_jacobian ( Elem & elem, uint side, const NumericVector<Number> & soln, 
                               SparseMatrix<Number> * jacobian , NumericVector<Number> * residual );

  virtual void set_bc( const RealTensor & sigtot_ ) { sigtot = sigtot_ ; }
  inline bool next_qp( bool inc=true );       /// QP control

  void init_fem();

  const MaterialConfig * config;

private:
  uint QP;                     /// The current quadrature point during integration
  optional<RealTensor> sigtot;
  ViscoplasticSolver & vpsolver; 
  TransientNonlinearImplicitSystem & system;

  /// FEM Stuff
  QGauss qrule;
  const Elem * elem;                               /// The element that the material has been reinit'ed to
  DenseVector<Number> Re;                          /// RHS vector for the element
  vector<dof_id_type> dof_indices;
  unique_ptr<FEBase> fe;                           /// The finite element object to hold shape funtions, jxw, etc
                       
  /// n_dofsv: number of dofs for each variable
  /// n_dofs: total number of dofs (all vars)
  /// n_uvars: Number of displacement vars. EG: U_x,y,z and Ueg_x,y,z
  uint n_dofs, n_dofsv, n_uvars, n_dofs_eg;

  string name;
  suint sid;   /// Subdomain id

private:
  /// Autodiff variables!
  AD::Vec ad_Uib, ad_Fib; // Flattened autodiff stuff
  inline void _init_autodiff() 
  { 
    uint ndofs = n_dofsv * 3 + n_dofs_eg*(n_uvars-3);
    ad_Uib.resize( ndofs ); ad_Uib.setZero();
    ad_Fib.resize( ndofs ); ad_Fib.setZero();
  }
  // 0,1,2 = Ux,Uy,Uz ; 3,4,5:UegX,UegY,UegZ
  inline uint Ndof(uint i) { return i>2 ?  n_dofs_eg : n_dofsv; }
  inline uint _idx( uint i, uint B ) { uint ret = i*n_dofsv + B; if ( i > 2 ) ret += (i-3) * n_dofs_eg; return ret; }
  inline AD::real & Uib( uint i, uint B ) { return ad_Uib[ _idx(i,B) ]; }
  inline AD::real & Fib( uint i, uint B ) { return ad_Fib[ _idx(i,B) ]; }
  //////

  AD::Vec residual_qp( const AD::Vec & /* ad_Uib */ );   
};

/**
 *  Advances in the quadrature integration point.
 *  Updates the Properties and temperature coupler pointer.
 */
inline bool ViscoPlasticMaterialBC::next_qp( bool inc )
{ 
  if ( inc ) 
    if ( ++QP == qrule.n_points() ) return false;
  return true; 
}

}} // ns
