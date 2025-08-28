
#pragma once

/**
 *
 *
 *
 */

#include "harpy/Global.h"
#include "config/InoutConfig.h"

namespace config {
namespace reader {

/* */
class InoutReader
{
public:
  InoutReader( InoutConfig & config );

  void check_files();
  void parse_inout_file();

  enum State { INITIAL, PROBE };

  InoutConfig & config;
  string filename;
  ProbeConfig * curr_probe_config;

  uint ln;               /// line number
  string line;           /// line being parsed
  State current_state;   /// Current state of the state machine
                         
  bool next_state();
  void probe_state();
};


}} // ns
