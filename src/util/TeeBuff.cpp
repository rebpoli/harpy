#include "util/TeeBuff.h"

namespace util {

// Construct a streambuf which tees output to both input
// streambufs.
teebuf::teebuf(std::streambuf * sb1, std::streambuf * sb2)
  : _sb1(sb1) , _sb2(sb2) 
{ }

// This tee buffer has no buffer. So every character "overflows"
// and can be put directly into the teed buffers.
int teebuf::overflow(int c)
{
  if ( RANK >= 1 ) return EOF; // só no proc 1
  if (c == EOF) { return !EOF; }
  else {
    if (!_sb1 || !_sb2) throw std::runtime_error("Null streambuf in teebuf");
    int r1 = _sb1->sputc(c);
    int r2 = _sb2->sputc(c);
    return r1 == EOF || r2 == EOF ? EOF : c;
  }
}

// Sync both teed buffers.
int teebuf::sync()
{
  if (!_sb1 || !_sb2) throw std::runtime_error("Null streambuf in teebuf");
  int r1 = _sb1->pubsync();
  int r2 = _sb2->pubsync();
  return r1 == 0 && r2 == 0 ? 0 : -1;
}

} // ns
