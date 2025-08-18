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
#include "util/Autodiff.h"

#include "util/CsvFile.h"

#include <boost/serialization/access.hpp>

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
  ViscoPlasticMaterial( suint sid_, const MaterialConfig * config,
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
//  virtual void residual_and_jacobian_qp_0 ();
  // Add the element contribution to the global jac and res
  void residual_and_jacobian (  Elem & elem,
                                const NumericVector<Number> & soln,
                                SparseMatrix<Number> * jacobian , 
                                NumericVector<Number> * residual );

  /// Updates the interface data at the probe points
  void update_probes();

  /// 
  void apply_strain_initialization_method();
  void update_initial_strain();

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
  vector<vector<dof_id_type>> dof_indices_var;
  vector<vector< DenseSubMatrix<Number> >> Ke_var;

  /// Autodiff variables!
  AD::Vec ad_Uib, ad_Fib; // Flattened autodiff stuff
  AD::Mat ad_Jijbm;
  inline void _init_autodiff() 
  { 
    ad_Uib.resize( n_dofsv * 3 );  ad_Uib.setZero();
    ad_Fib.resize( n_dofsv * 3 ); ad_Fib.setZero();
    ad_Jijbm.resize( n_dofsv*3, n_dofsv*3); ad_Jijbm.setZero();
  }
  inline AD::real & Uib( uint i, uint B ) { return ad_Uib[ i*n_dofsv + B ]; }
  inline AD::real & Fib( uint i, uint B ) { return ad_Fib[ i*n_dofsv + B ]; }
  inline AD::real Jijbm( uint i, uint j, uint B, uint M)
  { return ad_Jijbm(i*n_dofsv + B, j*n_dofsv + M); }
  // The residual function
  AD::Vec residual_qp( const AD::Vec & /* ad_Uib */ );
  //////
 
  /// Interface to hold the viscoplastic properties (in and out)
public:
  ViscoplasticIFC vp_ifc;
protected:
  VPProps * P;   

  /// Child boundary condition material. Owned by this object because it inherits all its properties
  ViscoPlasticMaterialBC * bc_material;

  /// The parent solver and system to assemble
  ViscoplasticSolver & vpsolver; 
  TransientNonlinearImplicitSystem & system;

  /// The stress system to export
  ExplicitSystem & stress_system;

  uint n_dofs, n_dofsv;

public:
  /// The stress engine
  StressPostProc stress_postproc;

  CsvFile dfile; // Debugging file
  int res_jac_k; // Debugging counter
private:
    /* Serialization routines - polymorphic serialization. */
    friend class boost::serialization::access;
    template<class Ar> void serialize(Ar& ar, const unsigned /*version*/) 
    { ar & vp_ifc; }
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

