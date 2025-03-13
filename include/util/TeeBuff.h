#include "base/Global.h"

#include <streambuf>

using namespace std;

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
