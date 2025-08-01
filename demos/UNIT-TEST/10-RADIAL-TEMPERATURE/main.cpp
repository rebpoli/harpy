
#include "util/GridFile.h"

#include <iomanip>

using namespace std;
/**
 *
 *  Testing ...
 *
 */
int main() 
{
  GridRadialFile grid("temperature.csv.gz");

  double time = 13000;
  double z = 5600;

  ilog << "Interpolated:";
  ilog << setw(10) << "time"
       << setw(10) << "r"
       << setw(10) << "z"
       << setw(10) << "T";

//  for ( double r = 0 ; r < 25 ; r += 0.5 )
  // r = -1 qr:51.8938 qz:5587.84  T:343 
  time = -1;
  double r = 51.8938;
  z = 5587.84;
  {
    double T = grid.at(time, r, z);
    ostringstream ss;
    ss << setw(10)  << time;
    ss << setw(10)  << r;
    ss << setw(10)  << z;
    ss << setw(10)  << T;
    ilog << ss.str();
  }
}



