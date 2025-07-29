#pragma once

#include "base/Global.h"
#include "config/ModelConfig.h"
#include "harpy/Material.h"
#include "harpy/ExplicitMaterial.h"

#include "libmesh/explicit_system.h"
namespace libMesh { class System; }

/**
 *
 *
 */
class StressPostProc : public ExplicitMaterial
{
public:
  StressPostProc( Material * refmat, ExplicitSystem & sys_ );
  virtual ~StressPostProc();

  virtual void init_fem();

  virtual string hello() { return "StressPostProc."; }
};

