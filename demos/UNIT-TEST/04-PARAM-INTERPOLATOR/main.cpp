
#include "util/SpatialDataReader.h"
#include "base/HarpyInit.h"

/**
 *
 *
 * Reads a map of properties and interpolates.
 * 
 * The properties are in file POR.
 * The queries are in file QUERY.
 *
 *
 */


using namespace std;

int main (int argc, char ** argv)
{
  HarpyInit init( argc, argv ); // Read config, debug etc

  // Read data file and query point list
  SpatialDataReader por( "POR" );
  SpatialDataReader query( "QUERY" );

  ilog << "READ DATA POINTS: " << por;
  ilog << "Inverse Distance Weighted Interpolation:";
  for (auto & pt : query.get_points() ) 
  {
    // Do the interpolation
    double val = por.value_at( pt.x, pt.y, pt.z );

    ilog << "Point (" << pt.x << "," << pt.y << "," << pt.z << "): " << val;
  }

  return 0;
}
