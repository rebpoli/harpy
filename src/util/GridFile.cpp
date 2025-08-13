
#include "util/GridFile.h"
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
  double t, r, z, temp;

  DataEntry(const string& csv_line) 
  {
    CSVLine tokens(csv_line);
    t      =   tokens.next<double>();
    r      =   tokens.next<double>();
    z      =   tokens.next<double>();
    temp   =   tokens.next<double>();
  }
};

/**
 *
 */
GridRadialFile::GridRadialFile(const string & filename) 
{
  ilog1 << "Reading grid file '" << filename << "' ...";
  GzStream file(filename);
  string line;
  file.getline(line); // skip header

  vector<DataEntry> entries;

  while (file.getline(line)) 
  {
    DataEntry e(line);
    entries.push_back(e);
    t.push_back(e.t);
    r.push_back(e.r);
    z.push_back(e.z);
  }

  auto remove_duplicates = [](vector<double>& v) 
  {
    sort(v.begin(), v.end());
    v.erase(unique(v.begin(), v.end()), v.end());
  };
  remove_duplicates(t);
  remove_duplicates(r);
  remove_duplicates(z);

  Nt = t.size();
  Nr = r.size();
  Nz = z.size();

  grid.resize(Nt * Nr * Nz, nan(""));

  for (const auto& e : entries) 
  {
    size_t it = find(t.begin(), t.end(), e.t) - t.begin();
    size_t ir = find(r.begin(), r.end(), e.r) - r.begin();
    size_t iz = find(z.begin(), z.end(), e.z) - z.begin();
    grid[it * Nr * Nz + ir * Nz + iz] = e.temp;
  }
}

///** 
// *   Trilinear interpolation
// */
//double GridRadialFile::at(double qt, double qr, double qz) const 
//{
//  if ( qr < min_radius ) qr = min_radius;

//  auto find_index = [](const vector<double>& vec, double val) {
//    auto it = upper_bound(vec.begin(), vec.end(), val);
//    size_t i = max<size_t>(1, it - vec.begin()) - 1;
//    return min(i, vec.size() - 2);
//  };

//  size_t it = find_index(t, qt),
//         ir = find_index(r, qr),
//         iz = find_index(z, qz);


//  double dt = (qt - t[it]) / (t[it+1] - t[it]),
//         dr = (qr - r[ir]) / (r[ir+1] - r[ir]),
//         dz = (qz - z[iz]) / (z[iz+1] - z[iz]);


//  auto at = [&](size_t i, size_t j, size_t k) {
//    return grid[index(i, j, k)];
//  };


//  auto lerp = [](double a, double b, double w) {
//    return a * (1 - w) + b * w;
//  };


//  double c00 = lerp(at(it,   ir,   iz),     at(it,   ir+1, iz),     dr),
//         c01 = lerp(at(it,   ir,   iz+1),   at(it,   ir+1, iz+1),   dr),
//         c10 = lerp(at(it+1, ir,   iz),     at(it+1, ir+1, iz),     dr),
//         c11 = lerp(at(it+1, ir,   iz+1),   at(it+1, ir+1, iz+1),   dr);


//  double c0 = lerp(c00, c01, dz);
//  double c1 = lerp(c10, c11, dz);

//  return lerp(c0, c1, dt);
//}


double GridRadialFile::at(double qt, double qr, double qz) const
{
  if ( qr < min_radius ) qr = min_radius;

  // --- time: choose the previous time index (or the first if qt < t[0]) ---
  auto prev_time_index = [](const std::vector<double>& vec, double val) {
    auto it = std::upper_bound(vec.begin(), vec.end(), val);
    if (it == vec.begin()) return size_t(0);          // before first -> use first
    return size_t(it - vec.begin() - 1);               // previous entry (ok even if val >= last)
  };
  const size_t it = prev_time_index(t, qt);

  // --- r,z: clamp query to [first, last] then find cell (lower index) ---
  auto clamp = [](double v, double lo, double hi) {
    return std::max(lo, std::min(v, hi));
  };
  qr = clamp(qr, r.front(), r.back());
  qz = clamp(qz, z.front(), z.back());

  auto find_index = [](const std::vector<double>& vec, double val) {
    auto it = std::upper_bound(vec.begin(), vec.end(), val);
    size_t i = std::max<size_t>(1, it - vec.begin()) - 1;
    return std::min(i, vec.size() - 2); // ensure i+1 valid
  };
  const size_t ir = find_index(r, qr);
  const size_t iz = find_index(z, qz);

  // --- local fractions in r and z (bilinear in space only) ---
  const double dr = (qr - r[ir]) / (r[ir+1] - r[ir]);
  const double dz = (qz - z[iz]) / (z[iz+1] - z[iz]);

  auto at_cell = [&](size_t i, size_t j, size_t k) {
    return grid[index(i, j, k)];
  };
  auto lerp = [](double a, double b, double w) { return a * (1.0 - w) + b * w; };

  // values at chosen time 'it', bilinear in (r,z)
  const double c00 = lerp(at_cell(it, ir,   iz),   at_cell(it, ir+1, iz),   dr);
  const double c01 = lerp(at_cell(it, ir,   iz+1), at_cell(it, ir+1, iz+1), dr);

  return lerp(c00, c01, dz);
}

/** 
 *   Trilinear interpolation - convert the radial symmetry in cartesian
 */
double GridRadialFile::at(double qt, double qx, double qy, double qz) const 
{
  // Calculate the radius we want to query in the radial geometry
  double qr = sqrt( qx*qx + qy*qy );
  double T = at( qt, qr, qz );
//  dlog(1) << "qt:" << qt << " qr:"<< qr <<" qz:" <<qz << "  T:" << T;
  return T;

}
