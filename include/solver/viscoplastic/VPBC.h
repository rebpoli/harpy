#pragma once

#include "harpy/Global.h"
#include "config/BCConfig.h"
#include <set>
#include <map>
#include <vector>
#include "solver/common/TemporalBC.h"
#include "libmesh/system.h"

/**
 *
 * Resolution of the boundary constrains with respect to a libmesh system.
 * Resolves variables, boundary ids etc.
 *
 * Translation of BCConfig => BC for a given time and system.
 *
 */

namespace config { class BCConfig; }
namespace timeloop { class Timestep; }


namespace solver {
namespace viscoplastic {

using namespace libMesh;
using config::BCConfig;
using timeloop::Timestep;
using solver::common::TemporalBC;

class BC 
{
  public:
    /** A class to index the BC map **/
    class ElemSide {
      public:
        ElemSide( uint e, uint s ) : eid(e), side(s) {}
        uint eid, side;
        bool operator<(const ElemSide & other) const;
    };
    /** Items to store the resolved boundary conditions **/
    class DirichletItem {
      public:
        DirichletItem( sint bid_, uint vid_, double val_, string vname_, string bname_ ) :
                              bid(bid_), vid(vid_), val(val_), vname(vname_), bname(bname_) {}
        sint bid; uint vid; double val; string vname, bname;
    };
    /**    **/
    class PenaltyItem {
      public:
        PenaltyItem( sint bid_, uint vid_, double pen_K_, double pen_val_, string pen_name_ ) :
                    bid(bid_), vid(vid_), pen_K(pen_K_), pen_val(pen_val_), pen_name(pen_name_) {}
        int bid; uint vid;  double pen_K, pen_val;  string pen_name;
    };
    /**    **/
    class ScalarItem {
      public:
        ScalarItem( sint bid_, uint vid_, uint svid_, string scalar_name_, string vname_, string bname_ ) :
                    bid(bid_), vid(vid_), svid(svid_), scalar_name(scalar_name_), vname(vname_), bname(bname_) {}
        sint bid; uint vid;  uint svid; string scalar_name, vname, bname;
    };
    /**    **/
    class STotItem {
      public:
        STotItem() : bid(0), val(0), bname("") {}
        STotItem( sint bid_, RealTensor val_, string bname_ ) :
                              bid(bid_), val(val_), bname(bname_) {}
        int bid;
        RealTensor val;
        string bname;
    };
    /** ** ** ** ** ** ** ** ** ** ** **/

    BC( const System & sys );
    ~BC() { _cleanup(); }
    bool update( double time );
    bool update( const Timestep & ts );

  private:
    const System & system;
    BCConfig & config;
    double time;        // Time of the current BC
    double reftime;    // Reference time of the current BC
    TemporalBC temporal_bcs;

    void _load_temporal_bcs();

    void _cleanup();
    void _validate();
    void _update_stot();
    void _update_dirichlet();
    void _update_scalar();
    void _update_penalty();
    void _update_temporal();

  public:
    vector< DirichletItem > dirichlet;
    vector< ScalarItem > scalar;

    map< ElemSide, STotItem *> stot;
    set< STotItem * > stotitem_ptrs; // List of created pointsr to manage cleanup

    map< ElemSide, vector<PenaltyItem *> > penalty;
    map< string, PenaltyItem * > penalty_ptrs; // List of created pointsr to manage cleanup

    friend Tester;
    friend ostream& operator<<(ostream& os, const BC & m);
};

ostream& operator<<(ostream& os, const BC & m);
ostream& operator<<(ostream& os, const BC::ElemSide & m);

ostream& operator<<(ostream& os, const vector<BC::DirichletItem> & m);
ostream& operator<<(ostream& os, const BC::DirichletItem & m);

ostream& operator<<(ostream& os, const vector<BC::ScalarItem> & m);
ostream& operator<<(ostream& os, const BC::ScalarItem & m);

ostream& operator<<(ostream& os, const map<BC::ElemSide, BC::STotItem *> & m);
ostream& operator<<(ostream& os, const BC::STotItem & m);

ostream& operator<<(ostream& os, const map<BC::ElemSide, vector<BC::PenaltyItem *>> & m);
ostream& operator<<(ostream& os, const vector<BC::PenaltyItem *> & m);
ostream& operator<<(ostream& os, const BC::PenaltyItem & m);

}} // ns
