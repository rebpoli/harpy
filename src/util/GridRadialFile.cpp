
#include "util/GridRadialFile.h"
#include "util/GzStream.h"

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>
#include <sstream>
#include <type_traits>

/**
 * Parse a CVS line and tokenize
 */
struct CSVLine 
{
  stringstream ss;
  CSVLine(const string& line) { ss.str(line); }

  template<typename T>
    T next() {
      string token;
      getline(ss, token, ',');
      if constexpr (is_same_v<T, int>) return stoi(token);
      else if constexpr (is_same_v<T, double>) return stod(token);
      else if constexpr (is_same_v<T, float>) return stof(token);
      else if constexpr (is_same_v<T, long>) return stol(token);
      else if constexpr (is_same_v<T, string>) return token;
      else flog << "Unsupported type.";
    }
};

/**
 *  Data entry for the interpolator
 */
struct DataEntry 
{
  double t, x, z, temp;

  DataEntry(const string& csv_line) 
  {
    CSVLine tokens(csv_line);
    t      =   tokens.next<double>();
    x      =   tokens.next<double>();
    z      =   tokens.next<double>();
    temp   =   tokens.next<double>();
  }
};

/**
 *
 */
GridRadialFile::GridRadialFile(const string & filename) 
{
  GzStream file(filename);
  string line;
  file.getline(line); // skip header

  vector<DataEntry> entries;

  while (file.getline(line)) 
  {
    DataEntry e(line);
    entries.push_back(e);
    t.push_back(e.t);
    x.push_back(e.x);
    z.push_back(e.z);
  }

  auto remove_duplicates = [](vector<double>& v) 
  {
    sort(v.begin(), v.end());
    v.erase(unique(v.begin(), v.end()), v.end());
  };
  remove_duplicates(t);
  remove_duplicates(x);
  remove_duplicates(z);

  Nt = t.size();
  Nx = x.size();
  Nz = z.size();

  grid.resize(Nt * Nx * Nz, nan(""));

  for (const auto& e : entries) 
  {
    size_t it = find(t.begin(), t.end(), e.t) - t.begin();
    size_t ix = find(x.begin(), x.end(), e.x) - x.begin();
    size_t iz = find(z.begin(), z.end(), e.z) - z.begin();
    grid[it * Nx * Nz + ix * Nz + iz] = e.temp;
  }
}

/** 
 *   Trilinear interpolation
 */
double GridRadialFile::at(double qt, double qx, double qz) const 
{
  auto find_index = [](const vector<double>& vec, double val) {
    auto it = upper_bound(vec.begin(), vec.end(), val);
    size_t i = max<size_t>(1, it - vec.begin()) - 1;
    return min(i, vec.size() - 2);
  };

  size_t it = find_index(t, qt),
         ix = find_index(x, qx),
         iz = find_index(z, qz);

  double dt = (qt - t[it]) / (t[it+1] - t[it]),
         dx = (qx - x[ix]) / (x[ix+1] - x[ix]),
         dz = (qz - z[iz]) / (z[iz+1] - z[iz]);

  auto at = [&](size_t i, size_t j, size_t k) {
    return grid[index(i, j, k)];
  };

  auto lerp = [](double a, double b, double w) {
    return a * (1 - w) + b * w;
  };

  double c00 = lerp(at(it,   ix,   iz),     at(it,   ix+1, iz),     dx),
         c01 = lerp(at(it,   ix,   iz+1),   at(it,   ix+1, iz+1),   dx),
         c10 = lerp(at(it+1, ix,   iz),     at(it+1, ix+1, iz),     dx),
         c11 = lerp(at(it+1, ix,   iz+1),   at(it+1, ix+1, iz+1),   dx);

  double c0 = lerp(c00, c01, dz);
  double c1 = lerp(c10, c11, dz);
  return lerp(c0, c1, dt);
}
