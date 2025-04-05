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

  // Poroelastic
  optional<double> porosity, permeability;
  optional<double> biot, bulk, skempton;

  // Files
  optional<string> porosity_file, permeability_file;
  optional<string> biot_file, bulk_file, skempton_file;

private:
  string line;
  uint ln;

  string model_dir ;
  string name, cfg;
  string filename;

  void check_files();
  void parse_material_file();

  void reg_param_str( string vname, string type, string val );
  void reg_param_dbl( string vname, string type, double val );
  optional<string> & _file_param( string & vname );
  optional<double> & _con_param( string & vname );

  friend ostream& operator<<(ostream& os, const MaterialConfig & m);
  bool operator<(const MaterialConfig & other) const;
};

ostream& operator<<(ostream& os, const MaterialConfig & m);
