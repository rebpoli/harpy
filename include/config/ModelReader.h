#pragma once

#include "base/Global.h"
#include "util/String.h"

#include <map>

/**
 *
 * Reads all the model raw information into C++ datastructures.
 *
 */

class ModelConfig;

class ModelReader
{
public:
  ModelReader( ModelConfig & config_ );

private:
  ModelConfig & config;


  // Parsing stuff
  enum class State { INITIAL, SYSTEM, TIMESTEP };
  uint ln;               /// line number
  string line;           /// line being parsed
  State current_state;   /// Current state of the state machine
  
  // File manip
  void check_files();

  // Model state machine
  void read_model();
  bool next_state();
  void timestep_state();
  void system_state();

  friend ModelConfig;
};
