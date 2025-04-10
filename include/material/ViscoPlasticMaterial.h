#pragma once

#include "base/Global.h"
#include "harpy/Coupler.h"
#include "harpy/Material.h"
#include "config/ModelConfig.h"
#include <optional>

/**
 *
 */

namespace libMesh { class System; }
class ViscoPlasticMaterialBC;

#include "libmesh/transient_system.h"

/**
 *
 *  THE MATERIAL CLASS
 *
 */
class ViscoPlasticMaterial : public Material
{
public:
  ViscoPlasticMaterial( suint sid_, const MaterialConfig & config, TransientNonlinearImplicitSystem & sys_ );
  virtual ~ViscoPlasticMaterial();
  virtual void init_fem();

  void reinit( const Elem & elem, uint side=255 );
  void reinit( const NumericVector<Number> & soln, const Coupler & coupler, const Elem & elem, uint side=255 );
  virtual void jacobian (const NumericVector<Number> & soln, SparseMatrix<Number> & jacobian );
  virtual void residual (const NumericVector<Number> & soln, NumericVector<Number> & residual );

  virtual Material * get_bc_material();

  virtual bool is_bc() { return false; }

protected:
  void setup_variables();

  vector<vector<dof_id_type>> dof_indices_var;
  vector<vector< DenseSubMatrix<Number> >> Ke_var;
//  vector< DenseSubVector<Number> > Re_var;

  vector<vector<Number>> Uib; /// first index is the node b, the second is the dimension i
  vector<vector<Number>> Fib;
 
  TransientNonlinearImplicitSystem & system;

  // PARAMETERS
  double implicit;
  inline double C_ijkl(uint qp, uint i, uint j, uint k, uint l);

  // Parameters from couplers (from config)
  vector<double> lame_mu, lame_lambda, alpha_d;
  // Variables from couplers (external solvers)
  vector< double > temperature;

  // THE CHILD MATERIAL
  ViscoPlasticMaterialBC * bc_material;
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
  ViscoPlasticMaterialBC( suint sid_, const MaterialConfig & config_, TransientNonlinearImplicitSystem & sys_ );

  virtual bool is_bc() { return true; }

  virtual void jacobian (const NumericVector<Number> & soln, SparseMatrix<Number> & jacobian );
  virtual void residual (const NumericVector<Number> & soln, NumericVector<Number> & residual );

  virtual void set_bc( const RealTensor & sigtot_ ) { sigtot = sigtot_ ; }

private:
  optional<RealTensor> sigtot;

};


/**
 *
 */
inline double ViscoPlasticMaterial::C_ijkl(uint qp, uint i, uint j, uint k, uint l) {
  const auto kd = Math::kronecker_delta;  // From Global

  return lame_lambda[qp] * ( kd(i,j)*kd(k,l) ) 
           + lame_mu[qp] * ( kd(i, k)*kd(j, l) + kd(i, l)*kd(j, k) );
}
