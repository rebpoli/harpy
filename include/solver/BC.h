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
        uint vid;
        double val;
    };
    /**    **/
    class DirichletItem {
      public:
        DirichletItem( uint bid_, uint vid_, double val_, string vname_, string bname_ ) :
                              bid(bid_), vid(vid_), val(val_), vname(vname_), bname(bname_) {}
        int bid;
        uint vid;
        double val;
        string vname, bname;
    };
    /** ** ** ** ** ** ** ** ** ** ** **/

    BC( const System & sys );
    void update( double time );

  private:
    const System & system;
    BCConfig config;
    double time;        // Time of the current BC
    double reftime;    // Reference time of the current BC

    // (eid,sid) => vector< (vid, double) >
    map< ElemSide , vector<Item> > BCMap;
    vector< DirichletItem > Dirichlet;

    void _validate();
    void _update_dbl();
    void _update_dirichlet();

    friend Tester;
    friend ostream& operator<<(ostream& os, const BC & m);
};

ostream& operator<<(ostream& os, const BC & m);
ostream& operator<<(ostream& os, const BC::ElemSide & m);
ostream& operator<<(ostream& os, const BC::Item & m);
ostream& operator<<(ostream& os, const BC::DirichletItem & m);
ostream& operator<<(ostream& os, const vector<BC::DirichletItem> & m);
