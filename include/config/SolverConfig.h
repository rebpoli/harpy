#pragma once

#include "base/Global.h"
#include "config/SolverConfig.h"

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

class SolverConfig
{
public:
  /**  Structs **/
  struct MatConfig {
    string name, cfg;
    MatConfig( string & m, string c ) : name(m), cfg(c) {}

    // To be indexable
    bool operator<(const MatConfig & other) const {
      if (name != other.name) return name < other.name;
      return cfg < other.cfg;
    }
  };
  struct MatConfigMap : public map<string, MatConfig> { }; // subdom => (mat, config)

  /* This is a complete object. All parameters must be set at all times (no std::optional here). */
  struct Numerical {
    Numerical() : ls_pc("bjacobi"), ls_solver("cg"), ls_atol(1e-10), ls_rtol(1e-6) { }
    string ls_pc;         /// preconditioner for the linear solver
    string ls_solver;     /// linear solver engine

    double ls_atol; /// linear solver absolute tolerance
    double ls_rtol; /// linear solver absolute tolerance
  };

  /** Data structure to the outside **/
  MatConfigMap material_config; /// the chosen configuration for this run
  Numerical numerical;

  SolverConfig( string model_dir_, string sys_name_, string sys_cfg_ );

  string model_dir, sys_file, sys_name, sys_cfg;

private:
  friend ostream& operator<<(ostream& os, const SolverConfig & m);
};

ostream& operator<<(ostream& os, const SolverConfig & m);
ostream& operator<<(ostream& os, const SolverConfig::MatConfigMap & m);
ostream& operator<<(ostream& os, const SolverConfig::MatConfig & m);
