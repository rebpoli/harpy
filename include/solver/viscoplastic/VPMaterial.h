#pragma once

#include "harpy/Global.h"


#include "config/ModelConfig.h"

#include "postproc/stress/StressPostProc.h"
#include "solver/viscoplastic/VPIfc.h"
#include "solver/thermal/TIfc.h"

#include "util/CsvFile.h"

#include "util/Autodiff.h"

#include <optional>
#include "libmesh/quadrature_gauss.h"
#include "libmesh/transient_system.h"

namespace solver {
namespace viscoplastic {

/**
 *
 */

class ViscoPlasticMaterialBC;
class ViscoplasticSolver;

using postproc::stress::StressPostProc;

/**
 *  THE MATERIAL CLASS
 */
class ViscoPlasticMaterial 
{
public:
  ViscoPlasticMaterial( suint sid_, const MaterialConfig * config,
                        TransientNonlinearImplicitSystem & sys_,
                        ViscoplasticSolver & vpsolver_,
                        bool called_from_bc_constructor = false);

  /**  Builds a material inheriting most properties from the reference.  */
//  ViscoPlasticMaterial( ViscoPlasticMaterial * refmat_ ) ;

  virtual ~ViscoPlasticMaterial();
  
  /// Initializes the FE objects. Assign quadrature to FEBase etc.
  void init_fem();

  /// Fetches properties from the configuration file
  void init_properties();

  ViscoPlasticMaterialBC * get_bc_material();

  /// Applies an element (and side?) to the FE objects.
  void reinit( const Elem & elem_, uint side=255 );
  void reinit( const NumericVector<Number> & soln, const Elem & elem_, uint side=255 );

  // Add res and jac to the element Ke and Re
  virtual void residual_and_jacobian_qp ();
//  virtual void residual_and_jacobian_qp_0 ();
  // Add the element contribution to the global jac and res
  void residual_and_jacobian (  Elem & elem,
                                const NumericVector<Number> & soln,
                                SparseMatrix<Number> * jacobian , 
                                NumericVector<Number> * residual );

  /// Updates the interface data at the probe points
  void update_probes();

  /// Calculate the stress_system stuff
  void project_stress( Elem & elem_ );
  void props_at( VPProps & props, const Point & pt, const Elem * elem );

  /// Rewinds the interfaces to the beginning of the timeste (for TS cutting)
  void rewind( Elem & elem_ );

  /// Differentiates ViscoPlasticMaterial and ViscoPlasticMaterialBC in runtime
  virtual bool is_bc() { return false; }

  /// Provide pointers to the outside to feed the material with properties
  ViscoplasticIFC & get_viscoplastic_interface() { return vp_ifc; }

  // Debugging
  string hello() { return "ViscoPlasticMaterial"; }

private:
  /// Creates the variables in the stress system
  void setup_variables();

protected:

  /// Advances the quadrature point during an integration. 
  /// For initialization purposes, inc can be false so that the first time we see QP=0
  inline bool next_qp( bool inc = true);
  uint QP;  /// The current quadrature point during integration
            
  /// Runtime containers 
  vector<vector< DenseSubMatrix<Number> >> Ke_var;


  /// Autodiff variables!
  AD::Vec ad_Uib, ad_Fib; // Flattened autodiff stuff
  AD::Mat ad_Jijbm;
  inline void _init_autodiff() 
  { 
    uint ndofs = n_dofsv * 3 + n_dofs_eg*(n_uvars-3);
    ad_Uib.resize( ndofs ); ad_Uib.setZero();
    ad_Fib.resize( ndofs ); ad_Fib.setZero();
    ad_Jijbm.resize( ndofs, ndofs ); ad_Jijbm.setZero();
  }
  // 0,1,2 = Ux,Uy,Uz ; 3,4,5:UegX,UegY,UegZ
  inline uint Ndof(uint i) { return i>2 ?  n_dofs_eg : n_dofsv; }
  inline uint _idx( uint i, uint B ) { uint ret = i*n_dofsv + B; if ( i > 2 ) ret += (i-3) * n_dofs_eg; return ret; }
  inline AD::real & Uib( uint i, uint B ) { return ad_Uib[ _idx(i,B) ]; }
  inline AD::real & Fib( uint i, uint B ) { return ad_Fib[ _idx(i,B) ]; }
  inline AD::real Jijbm( uint i, uint j, uint B, uint M)
  { return ad_Jijbm( _idx(i,B), _idx(j,M) ); }
  // The residual function
  AD::Vec residual_qp( const AD::Vec & /* ad_Uib */ );
  //////
 
  void _setup_fem();
  
public:
  const MaterialConfig * config;

  ViscoplasticIFC vp_ifc;
  VPProps * P;   

  /// Child boundary condition material. Owned by this object because it inherits all its properties
  ViscoPlasticMaterialBC * bc_material;

  /// The parent solver and system to assemble
  ViscoplasticSolver & vpsolver; 
  TransientNonlinearImplicitSystem & system;

  /// The stress system to export
  ExplicitSystem & stress_system;

  DenseMatrix<Number> Ke; /// Jacobian for the element
  DenseVector<Number> Re; /// RHS vector for the element

  vector<dof_id_type> dof_indices;

  unique_ptr<FEBase> fe;  /// The finite element object to hold shape funtions, jxw, etc
  QGauss qrule;
  const Elem * elem;   /// The element that the material has been reinit'ed to

  /// n_dofsv: number of dofs for each variable
  /// n_dofs: total number of dofs (all vars)
  /// n_uvars: Number of displacement vars. EG: U_x,y,z and Ueg_x,y,z
  uint n_dofs, n_dofsv, n_uvars, n_dofs_eg;

  string name;
  suint sid;   /// Subdomain id

  /// The stress engine
  StressPostProc stress_postproc;

//  CsvFile dfile; // Debugging file
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
  ViscoPlasticMaterialBC( suint sid_, const MaterialConfig * config_,
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
 *  Advances in the quadrature integration point.
 *  Updates the Properties and temperature coupler pointer.
 */
inline bool ViscoPlasticMaterial::next_qp( bool inc )
{ 
  if ( inc ) 
    if ( ++QP == qrule.n_points() ) return false;

  P = & ( vp_ifc.get(QP) );

  return true; 
}

}} // ns
