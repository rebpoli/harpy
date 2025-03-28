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
 * A system file can have a number of "configurations" so that a same
 * model can be P10, P50 and P90, for example.
 *
 */

class SystemConfig
{
public:
  /**  Structs **/
  struct MatConfig {
    string subdomain, mat, cfg;
    MatConfig( string & s, string & m, string c ) : subdomain(s), mat(m), cfg(c) {}
  };
  struct MatConfigMap : public map<string, MatConfig> { }; // subdom => (mat, config)

  /** Data structure to the outside **/
  MatConfigMap material_config; /// the chosen configuration for this run

  SystemConfig( string model_dir_, string sys_name_, string sys_cfg_ );

  string model_dir, sys_file, sys_name, sys_cfg;

private:
  friend ostream& operator<<(ostream& os, const SystemConfig & m);
};

ostream& operator<<(ostream& os, const SystemConfig & m);
ostream& operator<<(ostream& os, const SystemConfig::MatConfigMap & m);
ostream& operator<<(ostream& os, const SystemConfig::MatConfig & m);
