#pragma once

#include "harpy/Global.h"

namespace config { class BCConfig; }

namespace solver {
namespace common {

  using config::BCConfig;

/**
 */
class TemporalBC
{
public:
  TemporalBC( BCConfig & config ) ;

};


}} // ns
