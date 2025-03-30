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

class BCConfig 
{
  public:
    /* Sub structures */
    template <typename T>
    struct Item {
        Item() : bname("undef"), vname("undef"), value() {}
        string bname, vname;
        T value;
    };
    using ItemDbl = Item<double>;
    using ItemStr = Item<string>;
    /** **/
    struct ItemTensor {
        ItemTensor() : bname("undef"), vname("undef"), value( 3, vector<double>(3,0) ) {}
        string bname, vname;
        vector< vector<double> > value;
    };
    /** **/
    struct TimeEntry {
        TimeEntry() : drained(0), 
                      has_temperature(0), has_pressure(0),
                      temperature(0), pressure(0) {}

        // Bname => item
        map<string, vector<ItemDbl>> dbl_bcs;
        map<string, vector<ItemStr>> scalar_bcs;
        map<string, vector<ItemStr>> penalty_bcs;
        map<string, ItemTensor> stot_bcs;
      
        bool drained, has_temperature, has_pressure;
        double temperature, pressure;

        // Helpers to fill this struct
        void add_numerical_bc( string bname, string vname, double val );
        void add_scalar_bc( string bname, string vname, string scalar_name );
        void add_penalty_bc( string bname, string vname, string scalar_name );
    };
    /** **/
    struct PenaltyBC { 
      PenaltyBC() : K(0), value(0) { flog << "This function is not expected to be called."; }
      PenaltyBC( double K_, double value_ ) : K(K_), value(value_) {}
      double K; double value; 
    };
    /** ** ** ** ** ** ** ** **/

    BCConfig();

    void build_penalty() ;
    void build_scalars() ;
    void build_initial() ;
    void build_bcs() ;

    double get_reftime( double time ) ;

    string sys_name;

    map<string, PenaltyBC> penalty;  // penalty_name => entry
    set<string> scalars;             // list of scalars
    
    map< string, double > initial_by_vname;
    map< double, TimeEntry > entry_by_time;

    void all_bnames( set<string> & ret ) ;

    friend Tester;
    friend ostream& operator<<(ostream& os, const BCConfig & m);
};

ostream& operator<<(ostream& os, const BCConfig & m);
ostream& operator<<(ostream& os, const BCConfig::PenaltyBC & m);
ostream& operator<<(ostream& os, const map<string,BCConfig::PenaltyBC> & m);
ostream& operator<<(ostream& os, const BCConfig::TimeEntry & m);
ostream& operator<<(ostream& os, const map<double,BCConfig::TimeEntry> & m);
ostream& operator<<(ostream& os, const BCConfig::ItemTensor & m);
