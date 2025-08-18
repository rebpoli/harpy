
#include "restart/Header.h"
#include "restart/Util.h"
#include "restart/File.h"

namespace restart {

/**
 *
 */
void Header::write() 
{ 
  if ( ! f->is_root() ) flog << "Should write only on the root process! (rank=0).";
  _write( f->os, MAGIC );
  _write( f->os, N );
}

/**
 *
 */
void Header::read( )   
{ 
  string _m; _read( f->is, _m ); 
  if ( _m != MAGIC ) flog << "Wrong magic number in file. Something went wrong.";
  _read( f->is, N );
}

}
