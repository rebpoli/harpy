#pragma once

#include "base/Global.h"

namespace restart {

class File;

/**
 *   A few headers to help
 */
struct Header {
  File * f;

  /// Number of entries in 
  long unsigned int N; 
  /// For checking if we are in sync
  string MAGIC = "ABCD";

  void write();
  void read();

};

}
