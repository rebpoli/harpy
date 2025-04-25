// C++ includes
#include <iostream>

#ifndef AUTODIFF_EIGEN_FOUND 
#define AUTODIFF_EIGEN_FOUND ON
#endif

#include "util/Autodiff.h"

using namespace std;

// The vector function for which the Jacobian is needed
AD::Vec f(const AD::Vec & x)
{
  std::cout << endl << "f. x=" << endl << x << std::endl;
  AD::Mat A(5,5);

  A.row(1) = x * x.sum();
  return A.row(2); // x*x.sum(); // A(1,1);
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
