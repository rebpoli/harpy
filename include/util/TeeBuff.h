#include "harpy/Global.h"

#include <streambuf>

using namespace std;

namespace util {

class teebuf : public std::streambuf
{
  public:
    teebuf(std::streambuf * sb1, std::streambuf * sb2);
    int overflow(int c);
    int sync();

  private:
    std::streambuf* _sb1;
    std::streambuf* _sb2;
};

} // ns
