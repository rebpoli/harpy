#include "util/OutputOperators.h"

#include "libmesh/system.h"
#include "libmesh/point.h"

namespace util {

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

ostream& operator<<(ostream& os, const vector<optional<double>> & m)
{
  os << "[";
  for ( uint i=0; i<m.size() ; i++) {
    if ( i ) os << " ";
    os << m[i];
  }
  os << "] ("<< m.size() <<")";
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

ostream& operator<<(ostream& os, const util::CIMap<double> & m) {
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

}
