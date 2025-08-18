#pragma once

#include "base/Global.h"
#include <fstream>
#include <boost/mpi.hpp>

// Fwd declarations
class SLViscoplastic;
class ViscoplasticSolver;
class ViscoPlasticMaterial;
class Solver;

namespace restart {

class Header;

/**
 *
 *
 */
class File
{
  boost::mpi::communicator world;
  string filename;

  ofstream os;
  ifstream is;

  public:
    File( string filename_ ) : world(), filename(filename_) {};

    // Read/write of harpy entities
    void write( const SLViscoplastic & sloop );
    void write( const Solver * svr );
    void write( const ViscoplasticSolver * svr );

    void read( const SLViscoplastic & sloop );
    void read( const Solver * svr );
    void read( const ViscoplasticSolver * svr );

    inline bool is_root() { return world.rank() == 0; }
  private:

    friend Header;
};


}
