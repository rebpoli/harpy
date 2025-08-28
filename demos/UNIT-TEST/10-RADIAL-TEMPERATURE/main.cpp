
#include "util/GridFile.h"

#include <iomanip>

using namespace std;
using namespace util;

/**
 *
 *  Testing ...
 *
 */
int main() 
{
  GridRadialFile grid("temperature.csv.gz");

  double time = 400;
  double z = 5485;
  double r = 3;

  ilog << "Interpolated:";
  ilog << setw(10) << "time"
       << setw(10) << "r"
       << setw(10) << "z"
       << setw(10) << "T";

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



