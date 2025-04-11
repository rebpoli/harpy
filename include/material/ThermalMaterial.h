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

#include "libmesh/explicit_system.h"

/**
 *
 *  THE MATERIAL CLASS
 *
 */
class ThermalMaterial : public Material
{
public:
  ThermalMaterial( suint sid_, const MaterialConfig & config, ExplicitSystem & sys_ );
  virtual ~ThermalMaterial();
  virtual void init_fem();

  virtual void project( ElemCoupler & ec, string vname );
  void reinit( const Elem & elem, uint side=255 );
//  void reinit( const NumericVector<Number> & soln, const Coupler & coupler, const Elem & elem, uint side=255 );
//  virtual void jacobian (const NumericVector<Number> & soln, SparseMatrix<Number> & jacobian );
//  virtual void residual (const NumericVector<Number> & soln, NumericVector<Number> & residual );

  virtual bool is_bc() { return false; }

protected:
  void setup_variables();

  ExplicitSystem & system;
};
