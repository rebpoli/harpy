#pragma once

#include "base/Global.h"
#include "config/SolverConfig.h"
#include "config/MaterialConfig.h"
#include "config/BCConfig.h"

#include <map>

/**
 *
 * Stores all the model raw information into C++ datastructures.
 *
 * The reading and parsing is done by ModelReader, SolverReader etc.
 *
 */

class ModelConfig
{
public:
  /** Data structure to the outside **/
  map< string, double > timestep;

  string systemloop;
  map< string, string > system_cfgid;        /// sys_name -> config name
  map< string, SolverConfig > systems;       /// System configurations
  map< string, MaterialConfig > materials;   /// Material configurations
  BCConfig boundary_config;
  ModelConfig( string model_dir_ );

  string model_dir, model_file;

private:
  friend ostream& operator<<(ostream& os, const ModelConfig & m);
};

ostream& operator<<(ostream& os, const ModelConfig & m);
