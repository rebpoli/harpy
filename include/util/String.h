#pragma once

#include "base/Global.h"

#include <regex>

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

