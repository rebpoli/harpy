#pragma once

#include "harpy/Global.h"
#include "util/String.h"

#include <map>

namespace config {

class ModelConfig;

namespace reader {

/**
 *
 * Reads all the model raw information into C++ datastructures.
 *
 */


class ModelReader
{
public:
  ModelReader( ModelConfig & config_ );

private:
  ModelConfig & config;


  // Parsing stuff
  enum class State { 
    INITIAL, EXTERNAL, SYSTEMLOOP, TIMESTEP, PENALTY, SCALAR, TIME
  };

  uint ln;               /// line number
  string line;           /// line being parsed
  State current_state;   /// Current state of the state machine
  double current_time;   /// Current timestamp being parsed
  
  // File manip
  void check_files();

  // Model state machine
  void parse_model_file();
  bool next_state();
  void timestep_state();
  void system_state();
  void time_state();
  void external_state();
  void scalar_state();
  void penalty_state();

  friend ModelConfig;
};


}} // ns
