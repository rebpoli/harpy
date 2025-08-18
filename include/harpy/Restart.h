#pragma once

#include "base/Global.h"

#include <boost/mpi.hpp>

#include <type_traits>
#include <fstream>

class SLViscoplastic;
class ViscoplasticSolver;
class ViscoPlasticMaterial;
class Solver;

namespace harpy { namespace Restart {

class Header;

/**
 *
 *
 */
class File
{
  boost::mpi::communicator comm;
  string filename;

  ofstream os;
  ifstream is;

  public:
    File( string filename_ ) : comm(), filename(filename_) {};

    // Read/write of harpy entities
    void write( const SLViscoplastic & sloop );
    void write( const Solver * svr );
    void write( const ViscoplasticSolver * svr );
    void write( const ViscoPlasticMaterial * mat );

    void read( const SLViscoplastic & sloop );
    void read( const Solver * svr );
    void read( const ViscoplasticSolver * svr );
    void read( const ViscoPlasticMaterial * mat );

    inline bool is_root() { return comm.rank() == 0; }
  private:

    friend Header;
};


/**
 *
 */
template <typename T>
void _write(ofstream& os, const T& v) {
    static_assert(is_trivially_copyable<T>::value, "_write<T> requires trivially copyable T");
    os.write(reinterpret_cast<const char*>(&v), sizeof(T));
}

template <typename T>
void _read(ifstream& is, T& v) 
{
    static_assert(is_trivially_copyable<T>::value, "_read<T> requires trivially copyable T");
    is.read(reinterpret_cast<char*>(&v), sizeof(T));
}

/// Functions - implementation in cpp
void _write(ofstream& os, const string& s) ;
void _read(ifstream& is, string& s) ;

/**
 *   A few headers to help
 */
struct Header {
  File * f;

  uint version; uint N; string S; 
  const string MAGIC = "chbkjsd";

  ///
  void write() 
  { 
    if ( ! f->is_root() ) flog << "Should write only on the root process! (rank=0).";
    _write( f->os, MAGIC );
    _write( f->os, version );
    _write( f->os, N );
    _write( f->os, S ); 
  }

  ///
  void read( )   
  { 
    string _m; _read( f->is, _m ); 
    if ( _m != MAGIC ) flog << "Wrong magic number in file. Something went wrong.";

    _read( f->is, version );
    _read( f->is, N );
    _read( f->is, S );  
  }

};

}}
