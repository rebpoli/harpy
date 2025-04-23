#pragma once

#include "base/Global.h"

#include "harpy/Material.h"
#include "config/ModelConfig.h"

#include "postproc/StressPostProc.h"

#include "material/ViscoplasticIFC.h"
#include "material/ThermalIFC.h"

#include "libmesh/transient_system.h"
#include <optional>

// Autodiff stuff
#include <autodiff/forward/real.hpp>
#include <autodiff/forward/real/eigen.hpp>
namespace ad = autodiff;

/**
 *
 */

class ViscoPlasticMaterialBC;
class ViscoplasticSolver;

/**
 *  THE MATERIAL CLASS
 */
class ViscoPlasticMaterial : public Material
{
public:
  ViscoPlasticMaterial( suint sid_, const MaterialConfig & config,
                        TransientNonlinearImplicitSystem & sys_,
                        ViscoplasticSolver & vpsolver_,
                        bool called_from_bc_constructor = false);
  virtual ~ViscoPlasticMaterial();
  
  /// Initializes the FE objects. Assign quadrature to FEBase etc.
  virtual void init_fem();

  /// Fetches properties from the configuration file
  void init_properties();

  ViscoPlasticMaterialBC * get_bc_material();

  /// Applies an element (and side?) to the FE objects.
  void reinit( const Elem & elem_, uint side=255 );
  void reinit( const NumericVector<Number> & soln, const Elem & elem_, uint side=255 );

  // Add res and jac to the element Ke and Re
  virtual void residual_and_jacobian_qp ();
  // Add the element contribution to the global jac and res
  void residual_and_jacobian (  Elem & elem,
                                const NumericVector<Number> & soln,
                                SparseMatrix<Number> * jacobian , 
                                NumericVector<Number> * residual );

  /// Updates the interface Qp. Consider if this should be done inside the interface (?)
  void update_ifc_qp();

  /// Calculate the stress_system stuff
  void project_stress( Elem & elem_ );

  /// Differentiates ViscoPlasticMaterial and ViscoPlasticMaterialBC in runtime
  virtual bool is_bc() { return false; }

  /// Provide pointers to the outside to feed the material with properties
  ViscoplasticIFC & get_viscoplastic_interface() { return vp_ifc; }
  ThermalIFC & get_thermal_interface() { return th_ifc; }

  // Debugging
  string hello() { return "ViscoPlasticMaterial"; }

private:
  /// Calc Cijkl from lame constants
  inline double C_ijkl(uint i, uint j, uint k, uint l);

  /// Creates the variables in the stress system
  void setup_variables();

protected:

  /// Advances the quadrature point during an integration. 
  /// For initialization purposes, inc can be false so that the first time we see QP=0
  inline bool next_qp( bool inc = true);
  uint QP;  /// The current quadrature point during integration
            
  /// Runtime containers 
  vector<vector<dof_id_type>> dof_indices_var;
  vector<vector< DenseSubMatrix<Number> >> Ke_var;

  /// Autodiff variables!
  ad::VectorXreal _ad_Uib, _ad_grad_uij, _ad_Fib; // Flattened autodiff stuff
  Eigen::MatrixXd _ad_Jijbm;
  inline void _init_autodiff( uint nB ) 
  { 
    _ad_Uib.resize( nB * 3 );  _ad_Uib.setZero();
    _ad_Fib.resize( nB * 3 ); _ad_Fib.setZero();
    _ad_grad_uij.resize( 3 * 3 );  _ad_grad_uij.setZero();
    _ad_Jijbm.resize( nB*3, nB*3); _ad_Jijbm.setZero();
  }
  inline ad::real & Uib( uint i, uint B ) { return _ad_Uib[ 3*B + i ]; }
  inline ad::real & Fib( uint i, uint B ) { return _ad_Fib[ 3*B + i ]; }
  inline ad::real & grad_u( uint i, uint j ) { return _ad_grad_uij[ 3*i + j ]; }
  inline void set_Jijbm( uint i, uint j, uint B, uint M, double val) { _ad_Jijbm(3*B + i, 3*M +j) = val; }

  ad::VectorXreal residual_qp( const ad::VectorXreal & /* _ad_Uib */ );
 
  /// Interface to hold the viscoplastic properties (in and out)
  ViscoplasticIFC vp_ifc;
  ViscoplasticIFC::Props * P;   

  /// Interface to hold the thermal properties (in and out)
  ThermalIFC th_ifc;
  ThermalIFC::Props * T;   // The viscoplastic properties from the configuration

  /// Child boundary condition material. Owned by this object because it inherits all its properties
  ViscoPlasticMaterialBC * bc_material;

  /// The parent solver and system to assemble
  ViscoplasticSolver & vpsolver; 
  TransientNonlinearImplicitSystem & system;

  /// The stress system to export
  ExplicitSystem & stress_system;

  /// The stress engine
  StressPostProc stress_postproc;

};

/**
 *
 *  THE MATERIAL BOUNDARY CLASS:
 *  Inherits as much as possible from the material.
 *  The FEM structures are one dimension lower
 *
 */
class ViscoPlasticMaterialBC : public ViscoPlasticMaterial
{
public:
  ViscoPlasticMaterialBC( suint sid_, const MaterialConfig & config_,
                          TransientNonlinearImplicitSystem & sys_, 
                          ViscoplasticSolver & vpsolver_ );

  virtual bool is_bc() { return true; }

  void residual_and_jacobian_qp();
  void residual_and_jacobian ( Elem & elem, uint side, const NumericVector<Number> & soln, 
                               SparseMatrix<Number> * jacobian , NumericVector<Number> * residual );

  virtual void set_bc( const RealTensor & sigtot_ ) { sigtot = sigtot_ ; }

private:
  optional<RealTensor> sigtot;

};


/**
 *
 */
inline double ViscoPlasticMaterial::C_ijkl( uint i, uint j, uint k, uint l) {
  const auto kd = Math::kronecker_delta;  // From Global

  return P->lame_lambda * ( kd(i,j)*kd(k,l) ) 
           + P->lame_mu * ( kd(i, k)*kd(j, l) + kd(i, l)*kd(j, k) );
}

/**
 *  Advances in the quadrature integration point.
 *  Updates the Properties and temperature coupler pointer.
 */
inline bool ViscoPlasticMaterial::next_qp( bool inc )
{ 
  if ( inc ) 
    if ( ++QP == qrule.n_points() ) return false;

  P = & ( vp_ifc.get(QP) );
  T = & ( th_ifc.get(QP) );

  return true; 
}

