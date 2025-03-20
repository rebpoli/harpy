#pragma once

#include "base/Global.h"
#include "config/BCConfig.h"

#include <set>
#include <map>
#include <vector>

#include "libmesh/system.h"

/**
 *
 * Resolution of the boundary constrains with respect to a libmesh system.
 * Resolves variables, boundary ids etc.
 *
 * Translation of BCConfig => BC for a given time and system.
 *
 */

using namespace libMesh;

class BC {
  public:
    /** A class to index the BC map **/
    class ElemSide {
      public:
        ElemSide( uint e, uint s ) : eid(e), side(s) {}
        uint eid, side;
        bool operator<(const ElemSide & other) const;
    };
    /** A class to hold a single BC definition **/
    class Item {
      public:
        Item( uint vid_, double val_ ) : vid(vid_), val(val_) {}
        uint vid, val;
    };

  private:
    const System & system;
    BCConfig config;
    double time;        // Time of the current BC
    double reftime;    // Reference time of the current BC

    // (eid,sid) => vector< (vid, double) >
    map< ElemSide , vector<Item> > BCMap;

  public:
    BC( const System & sys );
    void update( double time );

    friend Tester;
    friend ostream& operator<<(ostream& os, const BC & m);
};

ostream& operator<<(ostream& os, const BC & m);
ostream& operator<<(ostream& os, const BC::ElemSide & m);
ostream& operator<<(ostream& os, const BC::Item & m);
