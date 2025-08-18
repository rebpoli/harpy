
#include "restart/File.h"
#include "restart/Header.h"
#include "restart/Util.h"

#include "solverloop/SLViscoplastic.h"

namespace mpi  = boost::mpi;

namespace restart { 

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
  world.barrier();
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
  // Localize the svr->material_by_sid to the root processor
  //
  // svr->material_by_sid < uint , Material * >
  //        mat->by_elem  < uint , vector<VPProps> > VPPropsByElemMap
  //
  // Binary organization:
  //     HEADER( material_by_sid )
  //        HEADER( by_elem )
  //          HEADER( vector<VPProps> )
  //             VPProps
  MaterialBySidMap global_map;
  const MaterialBySidMap & local_map = svr->material_by_sid;
  local_map.localize_to_one( global_map );

  if ( is_root() ) 
  {
    Header h{ this, svr->material_by_sid.size(), "MATERIAL_MAP" };
    h.write();
  }



  SCOPELOG(1);
  dlog(1) << "VISCOPLASTICSOLVER";
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


} // ns
