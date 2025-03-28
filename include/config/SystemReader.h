#pragma once

#include "base/Global.h"
#include "util/String.h"
#include "config/SystemConfig.h"

#include <map>

/**
 *
 * Reads all the model raw information into C++ datastructures.
 *
 */


class SystemReader
{
public:

  SystemReader( SystemConfig & config_ );

private:
  map< string, SystemConfig::MatConfigMap > all_mat_cfgs;  /// sys_cfgid => { subdom => (mat,  mat_cfgid) }

  SystemConfig & config;

  // Parsing stuff
  enum class State { INITIAL, CONFIG };
  uint ln;               /// line number
  string line;           /// line being parsed
  State current_state;   /// Current state of the state machine
  string curr_sys_cfg;   /// System configuration being read
  
  // File manip
  void check_files();

  // System state machine
  void parse_sys_file();
  bool next_state();
  void material_state();

  friend SystemConfig;
  friend ostream& operator<<(ostream& os, const SystemReader & m);
};

ostream& operator<<(ostream& os, const SystemReader & m);
