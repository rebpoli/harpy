
#pragma once

/*
 *     Sets a couple of aliases to make the code more readable.
 */

#include <autodiff/forward/real/eigen.hpp>

namespace AD
{
  using real = autodiff::real;
  using Mat = Eigen::MatrixXd;
  using Vec = autodiff::VectorXreal;
  using autodiff::jacobian;
}
