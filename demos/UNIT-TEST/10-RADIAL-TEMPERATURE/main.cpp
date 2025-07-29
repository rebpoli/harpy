
#include "util/GridFile.h"

/**
 *
 *  Testing ...
 *
 */
int main() 
{
  GridRadialFile grid("temperature.csv.gz");
  double T = grid.at(912*24*60*60, 66, 5585);
  ilog << "Interpolated T: " << T-273;
}



