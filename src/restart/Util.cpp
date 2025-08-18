
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

}
