#pragma once

#include "base/Global.h"
#include "config/SolverConfig.h"
#include "config/MaterialConfig.h"
#include "config/BCConfig.h"
#include "util/String.h"

#include <map>

/**
 *
 * Stores all the model raw information into C++ datastructures.
 *
 * The reading and parsing is done by ModelReader, SolverReader etc.
 *
 */
using namespace harpy_string;

class ModelConfig
{
public:
  ModelConfig( string model_dir_ );

  string systemloop;
  CIMap<string > system_cfgid;        /// sys_name -> config name
  CIMap< SolverConfig > solvers;       /// System configurations
  map< string, MaterialConfig > materials;   /// Material configurations
  map< string, double > timestep;
  BCConfig boundary_config;

  string model_dir, model_file;

// Getters with validation
  SolverConfig * solver_config( string name );

private:
  friend ostream& operator<<(ostream& os, const ModelConfig & m);
};

ostream& operator<<(ostream& os, const ModelConfig & m);
ostream& operator<<(ostream& os, const map<string,MaterialConfig> & m);

/** Global variable with the model configuration */
extern ModelConfig * MODEL;
