#pragma once

#include "base/Global.h"
#include "harpy/Material.h"
#include "config/ModelConfig.h"
#include "postproc/StressPostProc.h"
#include <optional>

/**
 *
 */

namespace libMesh { class System; }
class ViscoPlasticMaterialBC;
class ViscoplasticSolver;
class Timestep;

#include "libmesh/transient_system.h"

/**  Holds the data to feed viscoplastic properties into the material **/
struct ViscoPlasticIFC
{
  struct Props 
  {
    double lame_mu, lame_lambda, alpha_d, beta_e;
    double creep_carter_a, creep_carter_q, creep_carter_n;
    // Gradient of the displacement
    RealTensor grad_u;
    // Stresses
    RealTensor sigtot, sigeff, deviatoric;
    double von_mises;
    RealTensor plastic_strain;
//    double plastic_strain_rate;
  };

  // Transposed version of the properties
  struct PropsTranspose
  {
    PropsTranspose( vector<Props> * by_qp ) {
      for ( auto & p : *by_qp ) {
        sigtot.push_back( p.sigtot );
        sigeff.push_back( p.sigeff );
        von_mises.push_back( p.von_mises );
        deviatoric.push_back( p.deviatoric );
      }
    }
    vector<RealTensor> sigtot, sigeff, deviatoric;
    vector<double> von_mises;
  };

  // Data
  map< uint, vector<Props> > by_elem;
  vector<Props> * by_qp;  // By qp
  bool valid;
  // Helpers
  void reinit( uint eid, uint nqp ) { 
    by_qp = &( by_elem[eid] );
    if ( ! by_qp->size() ) { by_qp->resize(nqp); valid=0; }
  }
  Props & get( uint qp ) { return (*by_qp)[qp]; }
  uint size() { return by_qp->size(); }

};

/**  Holds the data to feed temperature into the material **/
struct ThermalIFC
{
  struct Props 
  { 
    double temperature; 
  };

  // Data
  map< uint, vector<Props> > by_elem;
  vector<Props> * by_qp;  // By qp

  // Helpers
  void reinit( uint eid, uint nqp ) {
    by_qp = &( by_elem[eid] );
    if ( ! by_qp->size() ) by_qp->resize(nqp); 
  }
  Props & get( uint qp ) { return (*by_qp)[qp]; }
  uint size() { return by_qp->size(); }
};





/**
 *  THE MATERIAL CLASS
 */
class ViscoPlasticMaterial : public Material
{
public:
  ViscoPlasticMaterial( suint sid_, const MaterialConfig & config,
                        TransientNonlinearImplicitSystem & sys_,
                        ViscoplasticSolver & vpsolver_ );

  virtual ~ViscoPlasticMaterial();
  virtual void init_fem();
  void init_properties();

  void reinit( const Elem & elem_, uint side=255 );
  void reinit( const NumericVector<Number> & soln, const Elem & elem, uint side=255 );
  ViscoPlasticMaterialBC * get_bc_material();

  void fetch_from_coupler();

  // Add the element contribution to the global jac and res
  void residual_and_jacobian (  Elem & elem,
                                const NumericVector<Number> & soln,
                                SparseMatrix<Number> * jacobian , 
                                NumericVector<Number> * residual );
  // Add res and jac to the element Ke and Re
  virtual void residual_and_jacobian_qp ();

  /// Updates the interface Qp. Consider if this should be done inside the interface (?)
  void update_ifc_qp();

  /// Calculate the stress_system stuff
  void project_stress();

  string hello() { return "ViscoPlasticMaterial"; }

  virtual bool is_bc() { return false; }

  // Provide pointers to the outside to feed the material with properties
  ViscoPlasticIFC & get_viscoplastic_interface() { return vp_ifc; }
  ThermalIFC & get_thermal_interface() { return th_ifc; }

private:

  inline double C_ijkl(uint i, uint j, uint k, uint l);   /// Material stuff

protected:
  // Advances the quadrature point during an integration. 
  // For initialization purposes, inc can be false so that the first time we see QP=0
  inline bool next_qp( bool inc = true);
  uint QP;  /// The current quadrature point during integration
            
  void setup_variables();
            
  vector<vector<dof_id_type>> dof_indices_var;
  vector<vector< DenseSubMatrix<Number> >> Ke_var;

  /// AD variables!
  vector<vector<Number>> Uib; /// first index is the dimention i, the second is the node B
  vector<vector<Number>> Fib;
 
  // Interfaces to hold the properties that can be fed by other solvers or configuration
  ViscoPlasticIFC vp_ifc;
  ViscoPlasticIFC::Props * P;   

  ThermalIFC th_ifc;
  ThermalIFC::Props * T;   // The viscoplastic properties from the configuration

  // THE CHILD MATERIAL
  ViscoPlasticMaterialBC * bc_material;

  // The parent solver
  ViscoplasticSolver & vpsolver; 
  TransientNonlinearImplicitSystem & system;

  // The stress engine
  ExplicitSystem & stress_system;
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

/** OUTPUT STREAMS **/
ostream& operator<<(ostream& os, const ViscoPlasticIFC & m);
