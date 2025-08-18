#pragma once

#include "base/Global.h"

namespace restart {

class File;

/**
 *   A few headers to help
 */
struct Header {
  File * f;

  uint version; uint N; string S; 
  const string MAGIC = "ABCD";

  void write();
  void read();

};

}
