// C++ includes
#include <iostream>

#define AUTODIFF_EIGEN_FOUND ON
#include "util/Autodiff.h"

using namespace std;

// The vector function for which the Jacobian is needed
AD::Vec f(const AD::Vec & x)
{
  std::cout << endl << "f. x=" << endl << x << std::endl;
  return x * x.sum();
}

namespace TEST 
{
  void run()
  {
    AD::Vec x(5);                           // the input vector x with 5 variables
      x << 1, 2, 3, 4, 5;                         // x = [1, 2, 3, 4, 5]

      AD::Vec F;                              // the output vector F = f(x) evaluated together with Jacobian matrix below

      AD::Mat J = AD::jacobian(f, wrt(x), at(x), F); // evaluate the output vector F and the Jacobian matrix dF/dx

      std::cout << "F = \n" << F << std::endl;    // print the evaluated output vector F
      std::cout << "J = \n" << J << std::endl;    // print the evaluated Jacobian matrix dF/dx
  }
}

int main()
{ TEST::run(); }
