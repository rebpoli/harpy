#pragma once

#include "base/Global.h"
#include "util/String.h"
#include "config/SolverConfig.h"

#include <map>

/**
 *
 * Reads all the model raw information into C++ datastructures.
 *
 */


class SolverReader
{
public:

  SolverReader( SolverConfig & config_ );

private:
  map< string, SolverConfig::MatConfigMap > all_mat_cfgs;  /// sys_cfgid => { subdom => (mat,  mat_cfgid) }

  SolverConfig & config;

  // Parsing stuff
  enum class State { INITIAL, EXTERNAL, CONFIG, MESH, NUMERICAL, FEM };
  uint ln;               /// line number
  string line;           /// line being parsed
  State current_state;   /// Current state of the state machine
  string curr_sys_cfg;   /// System configuration being read

  
  // File manip
  void check_files();

  // System state machine
  void parse_sys_file();
  bool next_state();

  // States
  void config_state();
  void external_state();
  void numerical_state();
  void fem_state();

  // Helper
  void set_grid_type( string str ) ;
  void set_origin( double x, double y, double z ) ;

  friend SolverConfig;
  friend ostream& operator<<(ostream& os, const SolverReader & m);
};

ostream& operator<<(ostream& os, const SolverReader & m);
