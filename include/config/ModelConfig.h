#pragma once

#include "base/Global.h"

/**
 *
 * Reads all the model raw information into C++ datastructures.
 *
 */

class ModelConfig
{
public:
  ModelConfig( string sys_name_ );

private:
  void check_files();
  void read_model();

  string sys_name;
  string dir, model_file;
  
  friend ostream& operator<<(ostream& os, const ModelConfig & m);
};

ostream& operator<<(ostream& os, const ModelConfig & m);
