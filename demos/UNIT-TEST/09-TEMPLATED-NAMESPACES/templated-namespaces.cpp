
#include "base/HarpyInit.h"

#ifndef AUTODIFF_EIGEN_FOUND 
#define AUTODIFF_EIGEN_FOUND ON
#endif

#include "util/Autodiff.h"
#include "libmesh/tensor_value.h"
#include "libmesh/vector_value.h"

// Autodiff
struct N1 {
  using Real = AD::real; 
  using Vec = AD::Vec;
  using Mat = AD::Mat;
};

// Primitive
struct N2 {
  using Real = double; 
  using Vec = libMesh::RealVectorValue;
  using Mat = libMesh::RealTensor;
};

using namespace std;

class TestClass
{
public:
  template < typename N >
  void foo() {
    using Mat = typename N::Mat;
    auto m = is_same_v<N, N1> ? Mat(3,3) : Mat();
    m(1,1) = 2;
    dlog(1) << "Matrix: " << m;
  };


};

// Example usage
int main(int argc, char* argv[]) {
  HarpyInit init( argc, argv );

  dlog(1) << "Hello world!";

  TestClass tc;
  tc.foo<N1>();
  tc.foo<N2>();
}
