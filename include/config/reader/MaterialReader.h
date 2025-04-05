#pragma once

/**
 *
 *
 *
 */

#include "base/Global.h"

class MaterialConfig;

class MaterialReader
{
public:
  MaterialReader( MaterialConfig & _config );

private:

  MaterialConfig & config;

  void check_files();
  void parse_material_file();

  // State machine
  enum class State { INITIAL, ENGINE, FEM };
  uint ln;               /// line number
  string line;           /// line being parsed
  State current_state;   /// Current state of the state machine
  string curr_sys_cfg;   /// System configuration being read
                         
  bool next_state();
  void fem_state();
  void engine_state();

  void reg_param_str( string vname, string type, string val );
  void reg_param_dbl( string vname, string type, double val );
};
