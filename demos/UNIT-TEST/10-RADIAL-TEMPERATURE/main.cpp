
#include "util/Grid3D.h"

/**
 *
 *  Testing ...
 *
 */
int main() 
{
  Grid3D grid("temperature.csv.gz");
  double T = grid.at(912*24*60*60, 66, 5585);
  ilog << "Interpolated T: " << T-273;
}



