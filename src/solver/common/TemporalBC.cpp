
#include "solver/common/TemporalBC.h"

#include "config/BCConfig.h"

#include "util/GzStream.h"
#include "util/CSVLine.h"

#include "util/OutputOperators.h"
#include "util/String.h"

namespace solver {
namespace common {

using util::CSVLine;

using util::operator<<;

/**
 *
 */
TemporalBC::TemporalBC( BCConfig & config ) 
{
  SCOPELOG(1);
  if ( config.temporal_file.empty() ) return;

  dlog(1) << "Reading temporal file: " << config.temporal_file;

  GzStream file(config.temporal_file);

  string line;
  file.getline(line); // skip header

  // Read as map for fast lookup
  map< pair<string,string>, TBCEntry > en_map;
  while ( file.getline(line) ) 
  {
    CSVLine tokens(line);
    auto time  = tokens.next<double>();
    auto bname    = tokens.next<string>();
    auto vname   = tokens.next<string>();
    auto value = tokens.next<double>();

    pair<string,string> key = make_pair(bname, vname);
    en_map[key].bname = bname;
    en_map[key].vname = vname;
    en_map[key].tt.add( time, value );
  }

  // Store as vector for sequential access
  for ( const auto& [key, value] : en_map ) 
    entries.push_back( value );
}

/**
 *
 */
ostream& operator<<(ostream& os, const TBCEntry & m)
{
  os  << "TBCEntry{ "
      << "bname: " << setw(10) << left << m.bname 
      << "vname: " << setw(10) << left << m.vname  
      << "tt : " << setw(10) << left << m.tt  
      << " }" << endl;
  
  return os;
}


}} // ns
