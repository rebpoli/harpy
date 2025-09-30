#pragma once

#include "harpy/Global.h"
#include <map>

namespace util 
{

/**
 *
 * Implements a table of 2 doubles that can interpolate in between.
 * 
 */

class TimeTable
{
private:
    std::map<double, double> table;

public:
    double interpolate(double time) const;

    /** inliners **/
    inline void   add(double time, double value) 
        { table[time] = value; }

    inline size_t size() const 
        { return table.size(); }

    inline void   clear() 
        { table.clear(); }

    inline bool   hasTime(double time) const 
        { return table.find(time) != table.end(); }

    inline bool   remove(double time) 
        { return table.erase(time) > 0; }

    inline double getValue(double time) const 
    {
        auto it = table.find(time);
        if (it == table.end()) flog << "Error: Time " << time << " not found in table";
        return it->second;
    }

    friend ostream& operator<<(ostream& os, const TimeTable & m);
};

ostream& operator<<(ostream& os, const TimeTable & m);

} // ns
