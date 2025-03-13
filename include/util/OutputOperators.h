#pragma once

#include "base/Global.h"

#include "libmesh/vector_value.h"

#include <vector>
#include <set>

namespace libMesh { class Point; }
ostream& operator<<(ostream& os, const vector<libMesh::Point> & m);
ostream& operator<<(ostream& os, const vector<libMesh::Gradient> & m);
ostream& operator<<(ostream& os, const vector<pair<uint,uint>> & m);
ostream& operator<<(ostream& os, const vector<double> & m);
ostream& operator<<(ostream& os, const set<double> & m);
ostream& operator<<(ostream& os, const vector<uint> & m);


ostream& operator<<(ostream& os, const map<string, double> & m);

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
