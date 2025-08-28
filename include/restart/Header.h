#pragma once

#include "harpy/Global.h"

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

protected:
  void write();
  void read();
};

/// Crates the header and writes
struct HeaderWrite : public Header {
  HeaderWrite( File * f_ , long unsigned int N, string MAGIC ) :
    Header( f_, N, MAGIC ) { write(); }
};

/// Creates the header and reads
struct HeaderRead : public Header {
  HeaderRead( File * f_ , string MAGIC ) :
    Header( f_, 0, MAGIC ) { read(); }
};

}
