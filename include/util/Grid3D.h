#pragma once

#include "base/Global.h"

#include <vector>


/**
 *   The grid interpolator
 */
class Grid3D 
{
  vector<double> t, x, z, grid;
  size_t Nt, Nx, Nz;
  size_t index(size_t it, size_t ix, size_t iz) const { return it * Nx * Nz + ix * Nz + iz; }

public:
  Grid3D(const string & filename) ;
  double at(double qt, double qx, double qz) const;
};

