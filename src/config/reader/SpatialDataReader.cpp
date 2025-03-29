
#include "config/reader/SpatialDataReader.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <regex>

// Pt nested class implementation
SpatialDataReader::Pt::Pt(
    double x_, double y_, double z_, double val_)
    : x(x_), y(y_), z(z_), val(val_) {}

// SpatialDataReader class implementation
SpatialDataReader::SpatialDataReader( string fn_ ) : fn(fn_) 
{
  readFromFile();
}

/**
 *
 *
 */
void SpatialDataReader::readFromFile() 
{
  ifstream file(fn);
  if (!file.is_open()) flog << "Unable to open file: " << fn;

  ilog << "Reding file '" << fn << "' ...";

  // Reads format "x y z value"
  string NUMBER_PATTERN = R"(([-+]?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?))";
  string pat = R"(^\s*)" + NUMBER_PATTERN + R"(\s+)" + NUMBER_PATTERN + R"(\s+)" +
        NUMBER_PATTERN + R"(\s+)" + NUMBER_PATTERN + R"(\s*(?:#.*)?$)";
  regex data_pattern = regex(pat);

  string line;
  size_t line_number = 0;
  while (getline(file, line)) 
  {
    line_number++;

    if (line.empty() || line.find_first_not_of(" \t") == string::npos) continue;

    smatch matches;
    if ( ! regex_match(line, matches, data_pattern)) 
      flog << "Invalid format on line " << line_number << ": " << line << endl;

    double x   = stod(matches[1]);
    double y   = stod(matches[2]);
    double z   = stod(matches[3]);
    double val = stod(matches[4]);

    // Feed the rtree
    Pt pt( x, y, z, val );
    RTreePoint rtreePoint(x, y, z);
    rtree.insert(make_pair(rtreePoint, pt));
    // For debugging only
    points.push_back( pt );
  }

  dlog(1) << "Read " << points.size() << " data points.";
}

/**
 *
 * n_nb: number of neighbors used for the interpolation.
 *
 */
double SpatialDataReader::value_at( double x, double y, double z, uint n_nb )
{
  vector<pair<RTreePoint, Pt>> nearestNeighbors;
  rtree.query(
      bgi::nearest(  RTreePoint(x, y, z), n_nb  ), 
      back_inserter( nearestNeighbors           )
   );

  double totalWeight = 0.0;
  double weightedSum = 0.0;
  Pt queryPoint( x, y, z );
  for (auto & neighbor : nearestNeighbors) 
  {
    Pt pt = neighbor.second;
    double distance = queryPoint.dist( pt );
    if (distance < 1e-10) return pt.val;  
    double weight = 1.0 / distance;
    totalWeight += weight;
    weightedSum += weight * pt.val;
  }

  return weightedSum / totalWeight;
}

/**
 *
 */
double SpatialDataReader::Pt::dist(const SpatialDataReader::Pt & m)
{
  return sqrt( pow(x - m.x, 2) + pow(y - m.y, 2) + pow(z - m.z, 2));
}

/**
 *
 *
 */
ostream& operator<<(ostream& os, const SpatialDataReader & m)
{
    os << "Total points read: " << m.points.size() << endl;
    for (auto pt : m.points) 
      os << pt;
    return os;
}

ostream& operator<<(ostream& os, const SpatialDataReader::Pt & m)
{
    os << "("  << setw(5) << m.x;
    os << ", " << setw(5) << m.y;
    os << ", " << setw(5) << m.z;
    os << "): " << m.val << endl;
    return os;
}

