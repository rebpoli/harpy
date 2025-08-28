#pragma once

#include "harpy/Global.h"

#include <type_traits>
#include <fstream>

#include "libmesh/tensor_value.h"
#include "libmesh/point.h"

namespace restart {

using namespace libMesh;

/**
 *
 */
template <typename T>
void _write(ofstream& os, const T& v) {
    static_assert(is_trivially_copyable<T>::value, "_write<T> requires trivially copyable T");
    os.write(reinterpret_cast<const char*>(&v), sizeof(T));
}

template <typename T>
void _read(ifstream& is, T& v) 
{
    static_assert(is_trivially_copyable<T>::value, "_read<T> requires trivially copyable T");
    is.read(reinterpret_cast<char*>(&v), sizeof(T));
}

/// Functions - implementation in cpp
void _write(ofstream& os, const string& s) ;
void _read(ifstream& is, string& s) ;

/// Libmesh stuff
void _write(ofstream& os, const RealTensor & A) ;
void _read(ifstream& is, RealTensor & A) ;

void _write(ofstream& os, const Point & A) ;
void _read(ifstream& is, Point & A) ;
}
