#pragma once

#include "base/Global.h"
#include "config/SolverConfig.h"
#include "config/MaterialConfig.h"
#include "config/BCConfig.h"
#include "config/TimestepConfig.h"
#include "util/String.h"

#include <map>

/**
 *
 * Stores all the model raw information into C++ datastructures.
 *
 * The reading and parsing is done by ModelReader, SolverReader etc.
 *
 */
using harpy_string::CIMap;

class ModelConfig
{
public:
  ModelConfig( string model_dir_ );

  string systemloop;
  CIMap<string> system_cfgid;        /// sys_name -> config name
  CIMap< SolverConfig > solvers;       /// System configurations
  set< MaterialConfig > materials;   /// Material configurations
  TimestepConfig timestep;
  BCConfig boundary_config;

  string model_dir, model_file;

// Getters with validation
  SolverConfig * solver_config( string name );

private:
  friend ostream& operator<<(ostream& os, const ModelConfig & m);
};

ostream& operator<<(ostream& os, const ModelConfig & m);
ostream& operator<<(ostream& os, const CIMap<MaterialConfig> & m);
ostream& operator<<(ostream& os, const TimestepConfig & m);

/** Global variable with the model configuration */
extern ModelConfig * MODEL;
