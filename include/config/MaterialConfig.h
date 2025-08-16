#pragma once

#include "base/Global.h"
#include <optional>
#include <map>
#include <vector>

#include "libmesh/point.h"
using namespace libMesh;

/**
 *
 * This class holds all the information for a given material.
 * Note that this is a physical description, so it may cover a number of systems.
 *
 * The material, however, does not need to be complete.
 * std::optional type is used to indicate whether a quantity was set or not.
 *
 */

/** Some support structs ... **/
struct CreepMD {
  struct SS { double sig0, n, q, stretch; };
  struct TR { double sig0, m, c, alpha_w; };
  vector<SS> ss;
  vector<TR> tr;
  double etr, etr_n;
  CreepMD() : etr(0), etr_n(0) {}
};

/** Initialization methods **/
enum class MatInitializeMethod { GRAVITY, HYDROSTATIC };
struct MatInitialize {
  using enum MatInitializeMethod;
  MatInitialize() : method( GRAVITY ) {}
  MatInitializeMethod method;
};

/**
 *
 */
class MaterialConfig
{
public:
  MaterialConfig( const string & model_dir_, const string & name_, const string & cfg_ );

  // Dummy constructor to serve as key for set search
  MaterialConfig( const string & name_, const string & cfg_ ) : name(name_), cfg(cfg_) {};

  void get_property( vector<double> & ret, string pname, const vector<Point> & xyz, string context ) const;
  double get_property( string pname, const Point & xyz, string context ) const;

  string engine;   /// Poroelastic, viscoplastic, porothermoelastic ...

  // Poroelastic
  optional<double> porosity, permeability;
  optional<double> biot, young, poisson, density;
  // Thermal
  optional<double> beta_e, beta_d, alpha_d;

  // Secondary variables (computed from the primary above)
  optional<double> lame_mu, lame_lambda, bulk_modulus;

  optional<CreepMD> creep_md;

  MatInitialize initialize;

  // Files
  // Poroelastic
  optional<string> porosity_file, permeability_file;
  optional<string> biot_file, young_file, poisson_file, density_file;
  // Thermal
  optional<string> beta_e_file, beta_d_file;

  // Get param reference by name
  optional<string> & file_param( string & vname, string context );
  optional<double> & con_param( string & vname, string context );
  const optional<string> & file_param( string & vname, string context ) const;
  const optional<double> & con_param( string & vname, string context ) const;

  bool operator<(const MaterialConfig & other) const;

  // Basic stuff
  string model_dir ;
  string name, cfg;
  string filename;


private:

  friend ostream& operator<<(ostream& os, const MaterialConfig & m);
};

ostream& operator<<(ostream& os, const MaterialConfig & m);
ostream& operator<<(ostream& os, const CreepMD & m);
