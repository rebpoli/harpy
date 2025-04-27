#pragma once

#include "base/Global.h"
#include "harpy/ExplicitMaterial.h"
#include "config/ModelConfig.h"

/**
 *
 */

#include "libmesh/explicit_system.h"
namespace libMesh { class System; }

class ViscoPlasticMaterial;

/**
 *
 *  THE MATERIAL CLASS
 *
 */
class ThermalPostProc : public ExplicitMaterial
{
public:
  ThermalPostProc( ViscoPlasticMaterial & ref_material, ExplicitSystem & sys_ );
  virtual ~ThermalPostProc() {};
  virtual void init_fem();
};
