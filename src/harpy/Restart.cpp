
#include "harpy/Restart.h"

#include "solverloop/SLViscoplastic.h"

namespace mpi  = boost::mpi;

namespace harpy { namespace Restart {

/**
 * 
 *    WRITERS
 *
 */

/**
 *   This is a collective task:
 *        1) gather info from all processors
 *        2) write only on the first processor
 *
 *   The os is only open in rank=0
 */
void File::write( const SLViscoplastic & sloop )
{
  if ( is_root() ) os.open( filename , ios::binary );

  write( sloop.viscoplastic );

  // Must handle the close
  if ( is_root() ) os.close(); 

  // MPI sync
  comm.barrier();
}

/**
 *
 */
void File::write( const Solver * svr )
{ flog << "Must cast the solver into the right class befor writing."; }

/**
 *
 */
void File::write( const ViscoplasticSolver * svr )
{
  // TODO: Collect stuff from all procs into the root (collective)

  if ( ! is_root() ) return; // only writes in the root. 

  SCOPELOG(1);
  dlog(1) << "VISCOPLASTICSOLVER";

  Header h{ this, 5, 100, "abcd"};
  h.write();

  for ( uint i=0 ; i<10 ; i++ )
  {
    _write(os, i);
    _write(os, to_string(i));
  }
}

/**
 *
 */
void File::write( const ViscoPlasticMaterial * mat )
{
}

/**
 * 
 *    READERS
 *
 */


/**
 *
 */
void File::read( const SLViscoplastic & sloop )
{
  is.open( filename , ios::binary );

  read( sloop.viscoplastic );

  // Must handle the close
  is.close();
}

/**
 *
 */
void File::read( const Solver * svr )
{
}

/**
 *
 */
void File::read( const ViscoplasticSolver * svr )
{
  Header header{ this };
  header.read();

  dlog(1) << "  HEADER ver(" << header.version << ") N(" << header.N << ") S(" << header.S << ")" ;
  for ( uint i=0 ; i<10 ; i++ )
  {
    uint I;
    _read(is, I);

    string S;
    _read(is, S);


    dlog(1) << "  (" << i << ") " << S ;
  }
}

/**
 *
 */
void File::read( const ViscoPlasticMaterial * mat )
{
}


/**
 *
 *     LOWER LEVEL READ/WRITE OF NON-TRIVIALLY-COPYABLE TYPES
 *
 */

/**
 *
 */
void _write(ofstream& os, const string& s) 
{
    uint n = static_cast<uint>(s.size());
    os.write(reinterpret_cast<const char*>(&n), sizeof(n));
    os.write(s.data(), static_cast<streamsize>(n));
}

/**
 *
 */
void _read(ifstream& is, string& s) 
{
    uint n = 0;
    is.read(reinterpret_cast<char*>(&n), sizeof(n));
    s.resize(static_cast<size_t>(n));
    is.read(&s[0], static_cast<streamsize>(n));
}

}} // ns
