#pragma once

#include "base/Global.h"
#include "config/ModelConfig.h"
#include "harpy/Material.h"

#include "libmesh/explicit_system.h"
namespace libMesh { class System; }

/**
 *
 *
 */
class StressMaterial : public Material
{
public:
  StressMaterial( suint sid_, const MaterialConfig & config, ExplicitSystem & sys_ );
  virtual ~StressMaterial();

  virtual void init_fem();

  virtual void feed_coupler( ElemCoupler & ec );

  // Add to the system so that we can visualize in paraview
  virtual void project( ElemCoupler & ec, string vname );
  void reinit( const Elem & elem, uint side=255 );

  inline double C_ijkl(uint qp, uint i, uint j, uint k, uint l);

protected:
  void setup_variables();

  ExplicitSystem & system;

  // The information needed to compute the stresses, in every qp
  vector<double> lame_mu, lame_lambda, alpha_d;
  vector<double> temperature;

  vector< RealVectorValue > U;   ///  Usage: U[qp](i)
  vector< RealTensor > GRAD_U; ///  Usage: GRAD_U[qp](i,j)
};

/**
 *
 */
inline double StressMaterial::C_ijkl(uint qp, uint i, uint j, uint k, uint l) {
  const auto kd = Math::kronecker_delta;  // From Global

  return lame_lambda[qp] * ( kd(i,j)*kd(k,l) ) 
           + lame_mu[qp] * ( kd(i, k)*kd(j, l) + kd(i, l)*kd(j, k) );
}
