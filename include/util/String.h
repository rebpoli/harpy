#pragma once

#include "base/Global.h"
#include <regex>
#include <set>
#include <map>

#include <boost/algorithm/string.hpp>
namespace harpy_string { 

  // Some dirty local aliases so we can rewrite or remove
  // of this context if we wish

  using boost::to_lower;
  using boost::to_lower_copy;

  using boost::to_upper;
  using boost::to_upper_copy;

  using boost::trim;
  using boost::trim_left;
  using boost::trim_right;

  using boost::join;
  
  using boost::iequals;       // Case insensitive comparison

  using boost::replace_all;

  /**
   *
   */
  inline string remove_comments_and_trim (string & str) 
  {
    return regex_replace(
                    regex_replace( str, regex(R"(#.*$)"), "" ),          
                    regex(R"(^\s+|\s+$)"), 
                    ""
              ); 
  }

  /** Define a new type: case insensitive set **/
  struct ci_cmp {
    bool operator() ( const string & a , const string & b ) const
    { return boost::ilexicographical_compare( a, b ); }
  };
  using CISet = set< string, ci_cmp >;
  template<typename T>
  using CIMap = map<string, T, ci_cmp>;

}
