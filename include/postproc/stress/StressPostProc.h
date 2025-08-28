#pragma once

#include "harpy/Global.h"
#include "config/ModelConfig.h"
#include "solver/common/ExplicitMaterial.h"

#include "libmesh/explicit_system.h"
namespace libMesh { class System; }

namespace postproc {
namespace stress {

using solver::common::ExplicitMaterial;
using solver::viscoplastic::ViscoPlasticMaterial;

using namespace libMesh;

/**
 *
 *
 */
class StressPostProc : public ExplicitMaterial
{
public:
  StressPostProc( ViscoPlasticMaterial * refmat, ExplicitSystem & sys_ );
  virtual ~StressPostProc();

  virtual void init_fem();

  virtual string hello() { return "StressPostProc."; }
};

}} // ns
