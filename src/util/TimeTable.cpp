#include "util/TimeTable.h"

namespace util {

/**
 *
 */
double TimeTable::interpolate(double time) const
{
  if (table.empty()) flog << "Error: TimeTable is empty";
  if (table.size() == 1) 
    return table.begin()->second;

  auto upper = table.lower_bound(time);
  if (upper != table.end())
  if (upper->first == time) 
    return upper->second;
  if (upper == table.begin()) 
    return table.begin()->second;

  if (upper == table.end()) 
    return table.rbegin()->second;

  auto lower = std::prev(upper);
  double t1 = lower->first, v1 = lower->second,
         t2 = upper->first, v2 = upper->second;
  return v1 + (v2 - v1) * (time - t1) / (t2 - t1);
}

/** **/
ostream& operator<<(ostream& os, const TimeTable & m)
{
  os << "TimeTable[" << m.table.size() << "]: ";
  for (const auto& pair : m.table)
    os << "(" << pair.first << "," << pair.second << ") ";
  return os;
}

} // ns
