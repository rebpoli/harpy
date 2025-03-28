#pragma once

#include "base/Global.h"
#include "config/SystemConfig.h"

#include <map>

/**
 *
 * Stores all the model raw information into C++ datastructures.
 *
 * The reading and parsing is done by ModelReader.
 *
 */

class ModelConfig
{
public:
  /** Data structure to the outside **/
  map< string, double > timestep;

  map< string, string > system_cfgid;   /// sys_name -> config name
  map< string, SystemConfig > systems;  /// Reads a system configuration

  ModelConfig( string model_dir_ );

  string model_dir, model_file;

private:
  friend ostream& operator<<(ostream& os, const ModelConfig & m);
};

ostream& operator<<(ostream& os, const ModelConfig & m);
