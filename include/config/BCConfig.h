#pragma once

#include "base/Global.h"

#include <set>
#include <map>
#include <vector>

// Numero de variaveis para indexacao
enum { SIGEFF=0, SIGTOT, HEAT };

/**
 *
 * Parses the configuration file. This class is unaware of the mesh.
 * It organizes the input data into datastructures.
 *
 */

class BCConfig {
  public:
    /* Subclasses */
    
    template <typename T>
    class Item {
      public:
        Item() : bname("undef"), vname("undef"), value() {}
        string bname, vname;
        T value;
    };
    using ItemDbl = Item<double>;
    using ItemStr = Item<string>;

    class ItemTensor {
      public:
        ItemTensor() : bname("undef"), vname("undef"), value( 3, vector<double>(3,0) ) {}
        string bname, vname;
        vector< vector<double> > value;
    };

    class TimeEntry {
      public:
        TimeEntry() : drained(0), 
                      has_temperature(0), has_pressure(0),
                      temperature(0), pressure(0) {}

        // Bname => item
        map<string, ItemDbl> dbl_bcs;
        map<string, ItemStr> str_bcs;
        map<string, ItemTensor> stot_bcs;
      
        bool drained, has_temperature, has_pressure;
        double temperature, pressure;
    };

    class PenaltyBC { public: double K; double value; };
    /* *** */

    BCConfig( string sys_name_ );

    void build_penalty() ;
    void build_scalars() ;
    void build_initial() ;
    void build_bcs() ;

    string sys_name;

    map<string, PenaltyBC> penalty;  // penalty_name => entry
    set<string> scalars;             // list of scalars
    
    map< string, double > initial_by_vname;
    map< double, TimeEntry > entry_by_time;

    friend Tester;
    friend ostream& operator<<(ostream& os, const BCConfig & m);
};

ostream& operator<<(ostream& os, const BCConfig & m);
ostream& operator<<(ostream& os, const BCConfig::PenaltyBC & m);
ostream& operator<<(ostream& os, const map<string,BCConfig::PenaltyBC> & m);
ostream& operator<<(ostream& os, const BCConfig::TimeEntry & m);
ostream& operator<<(ostream& os, const map<double,BCConfig::TimeEntry> & m);
ostream& operator<<(ostream& os, const BCConfig::ItemTensor & m);
