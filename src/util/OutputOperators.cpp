#include "util/OutputOperators.h"

#include "libmesh/system.h"
#include "libmesh/point.h"
ostream& operator<<(ostream& os, const vector<libMesh::Point> & m) {
  os << "[";
  for ( uint i=0; i<m.size() ; i++) {
    const libMesh::Point & p = m[i];
    if ( i ) os << " ";
    os << "(" << p(0) << "," << p(1) << "," << p(2) << ")";
  }
  os << "] ("<< m.size() <<")";
  return os;
}

ostream& operator<<(ostream& os, const vector<libMesh::Gradient> & m) {
  os << "[";
  for ( uint i=0; i<m.size() ; i++) {
    const libMesh::Point & p = m[i];
    if ( i ) os << " ";
    os << p;
  }
  os << "] ("<< m.size() <<")";
  return os;
}


ostream& operator<<(ostream& os, const vector<double> & m) {
  os << "[";
  for ( uint i=0; i<m.size() ; i++) {
    if ( i ) os << " ";
    os << m[i];
  }
  os << "] ("<< m.size() <<") {vec}";
  return os;
}

ostream& operator<<(ostream& os, const set<double> & m) {
  os << "[";
  uint i=0;
  for ( auto v : m ) {
    if ( i ) os << " ";
    os << v;
    i++;
  }
  os << "] ("<< m.size() <<") {set}";
  return os;
}

ostream& operator<<(ostream& os, const vector<pair<uint,uint>> & m)
{
  os << "[";
  for ( uint i=0; i<m.size() ; i++) {
    if ( i ) os << " ";
    os << m[i].first << "->" << m[i].second;
  }
  os << "] ("<< m.size() <<")";
  return os;
}

ostream& operator<<(ostream& os, const vector<uint> & m) {
  os << "[";
  for ( uint i=0; i<m.size() ; i++) {
    if ( i ) os << " ";
    os << m[i];
  }
  os << "] ("<< m.size() <<")";
  return os;
}

ostream& operator<<(ostream& os, const harpy_string::CIMap<double> & m) {
  os << "[";
  uint i=0;
  for ( auto v : m ) {
    if ( i++ ) os << ", ";
    os << "'" << v.first << "':";
    os << v.second;
  }
  os << "] ("<< m.size() <<")";
  return os;
}
ostream& operator<<(ostream& os, const map<string, double> & m) {
  os << "[";
  uint i=0;
  for ( auto v : m ) {
    if ( i++ ) os << ", ";
    os << "'" << v.first << "':";
    os << v.second;
  }
  os << "] ("<< m.size() <<")";
  return os;
}

ostream& operator<<(ostream& os, const pair<uint,uint> & m) {
  os << "("<<m.first<<","<<m.second<<")";
  return os;
}

ostream& operator<<(ostream& os, const map<pair<uint,uint>, bool> & m) {
  os << "[";
  uint i=0;
  for ( auto v : m ) {
    if ( i++ ) os << ", ";
    os << "'" << v.first << "':";
    os << v.second;
  }
  os << "] ("<< m.size() <<")";
  return os;
}

ostream& operator<<(ostream& os, const map<uint,uint> & m) {
  os << "[";
  uint i=0;
  for ( auto v : m ) {
    if ( i++ ) os << ", ";
    os << "'" << v.first << "':";
    os << v.second;
  }
  os << "] ("<< m.size() <<")";
  return os;
}
