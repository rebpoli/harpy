#pragma once

#include "base/Global.h"

#include <vector>

/**
 * 
 *  This is a virtual class providing a common interface.
 *
 */
class GridFile
{
public:
  GridFile() {}
  virtual ~GridFile() {}

  virtual double at(double t, double x, double y, double z) const 
  { flog << "Must be defined in the children"; return 0; };
};

/**
 *   The grid interpolator - radial symmetry
 */
class GridRadialFile : public GridFile
{
  vector<double> t, r, z, grid;
  size_t Nt, Nr, Nz;
  size_t index(size_t it, size_t ir, size_t iz) const { return it * Nr * Nz + ir * Nz + iz; }

public:
  GridRadialFile(const string & filename) ;

  double at(double qt, double qr, double qz) const;
  double at(double qt, double qx, double qy, double qz) const;
};

