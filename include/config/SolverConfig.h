#pragma once

#include "base/Global.h"
#include "base/Enums.h"
#include "config/SolverConfig.h"
#include "util/String.h"
#include "libmesh/point.h"

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

  /* FEM Stuff - variable order, family etc. */
  struct FEMSpec { 
    FEMSpec() : order("FIRST"), family("LAGRANGE"), type("CONTINUOUS") {};
    string order, family, type ;
    double implicit;
  };
  map<string,FEMSpec> fem_by_var;

  /* External: when the solution is known a priory and input as a file */
  struct ExternalFile {
    optional<string> filename;
    optional<GRID_TYPE> grid_type;
    optional<libMesh::Point> grid_origin;
    optional<double> min_radius;
    bool check();
    bool is_defined() { if (filename && grid_type && grid_origin) return true ; return false; }
  };

  /** Data structure to the outside **/
  MatConfigMap mat_config_by_name; /// the chosen configuration for this run
  Numerical numerical;

  SolverConfig( string model_dir_, string sys_name_, string sys_cfg_ );

  string model_dir, sys_file, sys_name, sys_cfg, mesh_filename;
  ExternalFile external_file;



private:
  friend ostream& operator<<(ostream& os, const SolverConfig & m);
};

ostream& operator<<(ostream& os, const SolverConfig & m);
ostream& operator<<(ostream& os, const SolverConfig::MatConfigMap & m);
ostream& operator<<(ostream& os, const SolverConfig::MatNameAndCfg & m);
ostream& operator<<(ostream& os, const SolverConfig::ExternalFile & m);
