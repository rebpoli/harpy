#pragma once

#include "base/Global.h"
#include <optional>

/**
 *
 * This class holds all the information for a given material.
 * Note that this is a physical description, so it may cover a number of systems.
 *
 * The material, however, does not need to be complete.
 * std::optional type is used to indicate whether a quantity was set or not.
 *
 */

class MaterialConfig
{
public:
  MaterialConfig( const string & model_dir_, const string & name_, const string & cfg_ );

  // Dummy constructor to serve as key for set search
  MaterialConfig( const string & name_, const string & cfg_ ) : name(name_), cfg(cfg_) {};

  string engine;   /// Poroelastic, viscoplastic, thermoporoelastic ...

  // Poroelastic
  optional<double> porosity, permeability;
  optional<double> biot, bulk, skempton;

  // Files
  optional<string> porosity_file, permeability_file;
  optional<string> biot_file, bulk_file, skempton_file;

  // FEM Stuff
  string fem_order, fem_family, fem_type;

  // Get param reference by name
  optional<string> & file_param( string & vname );
  optional<double> & con_param( string & vname );

  bool operator<(const MaterialConfig & other) const;

  // Basic stuff
  string model_dir ;
  string name, cfg;
  string filename;

private:

  friend ostream& operator<<(ostream& os, const MaterialConfig & m);
};

ostream& operator<<(ostream& os, const MaterialConfig & m);
