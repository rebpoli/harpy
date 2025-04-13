#pragma once

#include "base/Global.h"
#include "harpy/Material.h"
#include "config/ModelConfig.h"

/**
 *
 */

#include "libmesh/explicit_system.h"
namespace libMesh { class System; }

/**
 *
 *  THE MATERIAL CLASS
 *
 */
class ThermalMaterial : public MaterialExplicit
{
public:
  ThermalMaterial( suint sid_, const MaterialConfig & config, ExplicitSystem & sys_ );
  virtual ~ThermalMaterial();
  virtual void init_fem();

  void reinit( const Elem & elem, uint side=255 );

protected:
  void setup_variables();
};
