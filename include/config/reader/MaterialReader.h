#pragma once

/**
 *
 *
 *
 */

#include "harpy/Global.h"

namespace config {
class MaterialConfig;

namespace reader {


class MaterialReader
{
public:
  MaterialReader( MaterialConfig & _config );

private:

  MaterialConfig & config;

  void check_files();
  void parse_material_file();

  // State machine
  enum State { INITIAL, POROTHERMOELASTIC, CREEP, CREEP_MD, INITIALIZE };
  uint ln;               /// line number
  string line;           /// line being parsed
  State current_state;   /// Current state of the state machine
  string curr_sys_cfg;   /// System configuration being read
                         
  bool next_state();
  void porothermoelastic_state();
  void creep_md_state();
  void initialize_state();

  void reg_param_str( string vname, string type, string val, string context );
  void reg_param_dbl( string vname, string type, double val, string context );
};


}} // ns
