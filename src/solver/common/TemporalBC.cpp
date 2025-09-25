
#include "solver/common/TemporalBC.h"

#include "config/BCConfig.h"

#include "util/GzStream.h"
#include "util/CSVLine.h"

#include "util/OutputOperators.h"
#include "util/String.h"

namespace solver {
namespace common {

using util::CSVLine;

/**
 *  Data entry for the interpolator
 */
struct TBCEntry 
{
  double time, value;
  string bc, var;

  TBCEntry(const string& csv_line) 
  {
    CSVLine tokens(csv_line);
    time    =   tokens.next<double>();
    bc      =   tokens.next<string>();
    var     =   tokens.next<string>();
    value   =   tokens.next<double>();
  }

};
ostream& operator<<(ostream& os, const TBCEntry & m)
{
  os  << endl
      << "TBCEntry: { "
      << "time: " << setw(10) << left << m.time 
      << "bc: " << setw(10) << left << m.bc 
      << "var: " << setw(10) << left << m.var  
      << "value: " << setw(10) << left << m.value 
      << " }";
  
  return os;
}

using util::operator<<;

/**
 *
 *
 */
TemporalBC::TemporalBC( BCConfig & config ) 
{
  SCOPELOG(1);
  dlog(1) << "Reading temporal file: " << config.temporal_file;

  GzStream file(config.temporal_file);

  string line;
  file.getline(line); // skip header

  map< pair<string,string>, vector<TBCEntry> > entries;
  while ( file.getline(line) ) 
  {
    TBCEntry entry(line);
    pair<string,string> key = make_pair(entry.bc, entry.var);
    entries[key].push_back(entry);
  }

  dlog(1) << entries;

}


}} // ns
