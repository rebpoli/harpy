#pragma once

#include "base/Global.h"

#include "libmesh/vector_value.h"
#include "libmesh/tensor_value.h"
#include "libmesh/point.h"
#include "util/String.h"

#include <vector>
#include <set>
#include <optional>

namespace libMesh { class Point; }
ostream& operator<<(ostream& os, const vector<pair<uint,uint>> & m);
ostream& operator<<(ostream& os, const set<double> & m);


ostream& operator<<(ostream& os, const map<string, double> & m);
ostream& operator<<(ostream& os, const harpy_string::CIMap<double> & m);

ostream& operator<<(ostream& os, const pair<uint,uint> & m);
ostream& operator<<(ostream& os, const map<pair<uint,uint>, bool> & m);
ostream& operator<<(ostream& os, const map<uint,uint> & m);


/**
 * Rewriting the above little by little to take advantage of templates
 */
template <typename T>
ostream& operator<<(ostream& os, const set<T> & m)
{
  os << "[";
  uint i=0;
  for ( auto v : m ) {
    if ( i++ ) os << " ";
    os << "'" << v << "'";
  }
  os << "] ("<< m.size() <<")";
  return os;
}


/**
 *
 */
template <typename T>
ostream& operator<<(ostream& os, const vector<T> & m)
{
  os << "[";
  uint i=0;
  for ( auto v : m ) {
    if ( i++ ) os << " ";
    os << "'" << v << "'";
  }
  os << "] ("<< m.size() <<") {vec}";
  return os;
}


/*
 * A wrapper to override some operators provided by libraries to get a cleaner output
 *
 */
template<typename T>
struct Print {
  Print(const T& obj) : ref(obj) {}
  const T& ref;
};
template<typename T> ostream& operator<<(ostream& os, const Print<T>& printer) { os << printer.ref; return os; }

template<> 
inline ostream& operator<<(ostream& os, const Print<libMesh::Point>& printer)
{
  const libMesh::Point & p = printer.ref;
  os << "(";
  os << setw(3) << p(0) << ", ";
  os << setw(3) << p(1) << ", ";
  os << setw(3) << p(2) ;
  os << ")";
  return os;
}

template<> 
inline ostream& operator<<(ostream& os, const Print<vector<libMesh::Point>> & printer)
{
  os << "[";
  uint i=0;
  for ( auto v : printer.ref ) {
    if ( i++ ) os << " ";
    os << "'" << Print(v) << "'";
  }
  os << "] ("<< printer.ref.size() <<") {vec}";
  return os;
}

/* Output of optional stuff */
template <typename T>
ostream& operator<<(ostream& os, const optional<T> & m)
{
  if ( m ) 
    os << Print(*m) ; 
  else 
    os << "(undef)";
  return os;
}

