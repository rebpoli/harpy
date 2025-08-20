
#include "restart/Util.h"


namespace restart {

/**
 *
 */
void _write(ofstream& os, const string& s) 
{
    uint n = static_cast<uint>(s.size());
    os.write(reinterpret_cast<const char*>(&n), sizeof(n));
    os.write(s.data(), static_cast<streamsize>(n));
}

/**
 *
 */
void _read(ifstream& is, string& s) 
{
    uint n = 0;
    is.read(reinterpret_cast<char*>(&n), sizeof(n));
    s.resize(static_cast<size_t>(n));
    is.read(&s[0], static_cast<streamsize>(n));
}

// Write a 3x3 RealTensor
void _write(ofstream& os, const RealTensor& A)
{
    array<double, 9> buf;
    size_t k = 0;
    for (unsigned i = 0; i < 3; ++i)
        for (unsigned j = 0; j < 3; ++j)
            buf[k++] = A(i, j);

    os.write(reinterpret_cast<const char*>(buf.data()),
             static_cast<streamsize>(buf.size() * sizeof(buf[0])));
}

// Read a 3x3 RealTensor
void _read(ifstream& is, RealTensor& A)
{
    array<double, 9> buf{};
    is.read(reinterpret_cast<char*>(buf.data()),
            static_cast<streamsize>(buf.size() * sizeof(buf[0])));

    size_t k = 0;
    for (unsigned i = 0; i < 3; ++i)
        for (unsigned j = 0; j < 3; ++j)
            A(i, j) = buf[k++];
}

// Point
void _write(ofstream& os, const Point & A)
{
  _write( os, A(0) );
  _write( os, A(1) );
  _write( os, A(2) );
}

// Point
void _read(ifstream& is, Point & A)
{
  double x,y,z;
  _read( is, x );
  _read( is, y );
  _read( is, z );
  A = Point(x,y,z);
}


}
