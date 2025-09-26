#pragma once

#include "harpy/Global.h"
#include "util/TimeTable.h"
#include <vector>

/*
 *
 *      Provides access functions to interpolate time dependent boundary
 *      conditions.
 *
 *      The datastructure is a vector of entries, each describing a 
 *      Boundary (XP, XN, WELL, etc), a variable (UX, UY, SXX ...) 
 *      and a time table ( time: value ). 
 *
 */


namespace config { class BCConfig; }

namespace solver {
namespace common {

using config::BCConfig;
using util::TimeTable;

/** Datastructure */
struct TBCEntry 
{ 
  string bname, vname; // boundary and variable names
  TimeTable tt; 
};
ostream& operator<<(ostream& os, const TBCEntry & m);

/**
 */
class TemporalBC
{
public:
  TemporalBC( BCConfig & config ) ;

  vector<TBCEntry> entries;
};


}} // ns
