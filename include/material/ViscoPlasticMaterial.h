#pragma once

#include "base/Global.h"
#include "harpy/Material.h"
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
  ViscoPlasticMaterial( suint sid_, const MaterialConfig & config, System & sys_ );
  virtual ~ViscoPlasticMaterial();
  virtual void init_fem();

  void reinit( const NumericVector<Number> & soln, const Elem & elem, uint side=255 );
  virtual void jacobian (const NumericVector<Number> & soln, SparseMatrix<Number> & jacobian );
  virtual void residual (const NumericVector<Number> & soln, NumericVector<Number> & residual );

  virtual Material * get_bc_material();

  virtual bool is_bc() { return false; }

protected:
  void setup_variables();

  vector<vector<dof_id_type>> dof_indices_var;
  vector<vector< DenseSubMatrix<Number> >> Ke_var;
  vector< DenseSubVector<Number> > Re_var;

  vector<vector<Number>> Uib; /// first index is the node b, the second is the dimension i
  vector<vector<Number>> Fib;
 
  TransientNonlinearImplicitSystem & system;

  // PARAMETERS
  double implicit;
  inline double C_ijkl(uint i, uint j, uint k, uint l);

  // These material parameters are constant for now.
  // They must be imported from a Coupler.
  double lame_mu, lame_lambda; /// Lame parameters

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
  ViscoPlasticMaterialBC( suint sid_, const MaterialConfig & config_, System & sys_ );

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
inline double ViscoPlasticMaterial::C_ijkl(uint i, uint j, uint k, uint l) {
  const auto kd = Math::kronecker_delta;  // From Global

  return lame_lambda * ( kd(i,j)*kd(k,l) ) 
           + lame_mu * ( kd(i, k)*kd(j, l) + kd(i, l)*kd(j, k) );
}
