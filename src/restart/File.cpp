
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
void File::write( const ViscoplasticSolver * svr )
{
  SCOPELOG(1);
  // Localize the svr->material_by_sid to the root processor
  //
  // svr->material_by_sid < uint , Material * >
  //        mat->by_elem  < uint , vector<VPProps> > VPPropsByElemMap
  //
  // Binary organization:
  //     HEADER( material_by_sid )       -- in all procs; sid=subdomain_id
  //        vp_mat -> vp_ifc 
  //        HEADER( by_elem )            -- need to localize
  //          HEADER( vector<VPProps> )

  const map< uint, Material *> & material_by_sid = svr->material_by_sid;

  // Wrtie a header for teh map
  if ( is_root() ) 
    HeaderWrite( this, material_by_sid.size(), "MAT_BY_SID" );   /* WRITE */

  for ( auto & [ sid, mat ] : material_by_sid )
  {
    ViscoPlasticMaterial * vp_mat = dynamic_cast<ViscoPlasticMaterial *>(mat);
    ViscoplasticIFC & vp_ifc = vp_mat->vp_ifc;
    VPPropsByElemMap local_map = vp_ifc.by_elem;
    VPPropsByElemMap global_map;
    local_map.localize_to_one( global_map );

    // Now we have all data in rank=0 for this _sid_ . Write the file.
    if ( is_root() ) 
    {
      _write( os , sid );                                 /* WRITE */

      HeaderWrite(this, global_map.size(), "VEC_BY_EID"); /* WRITE */

      for ( auto & [eid, vec] : global_map ) 
      {
        _write( os, eid );                                /* WRITE */

        HeaderWrite(this, vec.size(), "EID_VEC");         /* WRITE */

        // Write all we want to store
        for ( VPProps & p : vec ) {
          _write( os, p.initial_strain );                 /* WRITE */
          _write( os, p.sigtot );                         /* WRITE */
        }
      }
    }
  }
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
void File::read( const ViscoplasticSolver * svr )
{
  SCOPELOG(1);
  HeaderRead h0( this, "MAT_BY_SID" );

  // Our target map to feed
  // NOTE: the map is build. We only feed selected data
  const map< uint, Material *> & material_by_sid = svr->material_by_sid;
  const MeshBase & mesh = svr->get_mesh();

  // SID LOOP
  for ( uint i=0 ; i< h0.N ; i++ )
  {
    uint sid; _read( is, sid );                                            /* READ */
    HeaderRead h1( this, "VEC_BY_EID" );                                   /* READ */

    // Feed the data read from file - only in the processor owning the element
    ViscoPlasticMaterial * vp_mat = dynamic_cast<ViscoPlasticMaterial *>( material_by_sid.at(sid) );
    ViscoplasticIFC & vp_ifc = vp_mat->vp_ifc;
    VPPropsByElemMap & by_elem = vp_ifc.by_elem;

    // EID LOOP
    for ( uint j=0; j<h1.N; j++ )
    {
      uint eid; _read( is, eid );                                          /* READ */
      HeaderRead h2( this, "EID_VEC" );                                    /* READ */

      bool is_local = 0;
      if ( mesh.processor_id() == mesh.elem_ref(eid).processor_id() ) is_local = 1;

      // Consistency check 
      if ( is_local ) {
        vector<VPProps> & vec = by_elem.at(eid);
        if ( vec.size() != h2.N ) flog << "Wrong vec size? " << vec.size() << " != " << h2.N;
      }

      // QUADRATURE POINT LOOP
      for ( uint qp=0; qp<h2.N; qp++ )
      {
        RealTensor initial_strain; _read(is, initial_strain);              /* READ */
        RealTensor sigtot;         _read(is, sigtot);                      /* READ */

        if ( ! is_local ) continue;  // it is important to do the dummy runs above to keep sync
                                     
        // Update the local structures
        vector<VPProps> & vec = by_elem.at(eid);
        VPProps &p = vec[qp];

        p.initial_strain = initial_strain;
      }
    }

  }
}



/** **/
void File::write( const Solver * svr )
{ flog << "Must cast the solver into the right class befor writing."; }
void File::read( const Solver * svr )
{ flog << "Must cast the solver into the right class befor writing."; }


} // ns
