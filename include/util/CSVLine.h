
#pragma once

#include "harpy/Global.h"

namespace util {

/**
 * Parse a CVS line and tokenize
 */
struct CSVLine 
{
  stringstream ss;
  CSVLine(const string& line) { ss.str(line); }

  template<typename T>
    T next() {
      string token;
      getline(ss, token, ',');
      if constexpr (is_same_v<T, int>) return stoi(token);
      else if constexpr (is_same_v<T, double>) return stod(token);
      else if constexpr (is_same_v<T, float>) return stof(token);
      else if constexpr (is_same_v<T, long>) return stol(token);
      else if constexpr (is_same_v<T, string>) return token;
      else flog << "Unsupported type.";
    }
};


} //ns
