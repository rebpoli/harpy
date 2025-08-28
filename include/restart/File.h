#pragma once

#include "harpy/Global.h"
#include <fstream>
#include <boost/mpi.hpp>

// Fwd declarations
namespace solverloop { class SLViscoplastic; } 
namespace solver { 
  namespace viscoplastic { class ViscoplasticSolver; class ViscoPlasticMaterial; } 
  namespace common { class Solver; }
}

namespace restart {

using solver::common::Solver;
using solver::viscoplastic::ViscoplasticSolver;
using solverloop::SLViscoplastic;

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

    void read( SLViscoplastic & sloop );
    void read( Solver * svr );
    void read( ViscoplasticSolver * svr );

    inline bool is_root() { return world.rank() == 0; }
  private:

    friend Header;
};


}
