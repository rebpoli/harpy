#pragma once

#include "base/Global.h"
#include "config/SolverConfig.h"
#include "util/String.h"

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

using harpy_string::CIMap;

class SolverConfig
{
public:
  /**  Structs **/
  struct MatNameAndCfg {
    string name, cfg;
    MatNameAndCfg( string & m, string c ) : name(m), cfg(c) {}

    // To be indexable
    bool operator<(const MatNameAndCfg & other) const {
      if (name != other.name) return name < other.name;
      return cfg < other.cfg;
    }
  };
  // subdom_name => (mat, config)
  struct MatConfigMap : public CIMap<MatNameAndCfg> { }; 

  /* This is a complete object. All parameters must be set at all times (no std::optional here). */
  struct Numerical {
    Numerical() : ls_pc("bjacobi"), ls_solver("cg"), ls_atol(1e-10), ls_rtol(1e-6) { }
    string ls_pc;         /// preconditioner for the linear solver
    string ls_solver;     /// linear solver engine

    double ls_atol; /// linear solver absolute tolerance
    double ls_rtol; /// linear solver absolute tolerance
  };

  /** Data structure to the outside **/
  MatConfigMap mat_config_by_name; /// the chosen configuration for this run
  Numerical numerical;

  SolverConfig( string model_dir_, string sys_name_, string sys_cfg_ );

  string model_dir, sys_file, sys_name, sys_cfg, mesh_filename;

private:
  friend ostream& operator<<(ostream& os, const SolverConfig & m);
};

ostream& operator<<(ostream& os, const SolverConfig & m);
ostream& operator<<(ostream& os, const SolverConfig::MatConfigMap & m);
ostream& operator<<(ostream& os, const SolverConfig::MatNameAndCfg & m);
