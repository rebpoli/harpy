
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

//     Now we have all data in rank=0 for this _sid_ . Write the file.
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
          _write( os, p.sigtot );                         /* WRITE */
        }
      }
    }

    /** Save probes **/
    //  HEADER ( probes_by_pname_by_elem )
    //    HEADER( by_pname )
    //      HEADER( by_elem )
    ProbeByPnameByElemMap & probes_by_pname_by_elem = vp_ifc.probes_by_pname_by_elem;
    ProbeByPnameByElemMap gl_map;
    vp_ifc.probes_by_pname_by_elem.localize_to_one( gl_map );

    if ( is_root() )
    {
      HeaderWrite( this, gl_map.size(), "PROBES_BYBY" );  /* WRITE */
      for ( auto &[ pname, by_elem ] : gl_map )
      {
        _write( os, pname );                                        /* WRITE */
        HeaderWrite( this, by_elem.size(), "PROBES_BY_ELEM" );      /* WRITE */
        for ( auto & [eid, vec] : by_elem ) 
        {
          _write( os, eid );                                        /* WRITE */
          dlog(1) << "WRITING EID:" << eid << " vec.size:" << vec.size();
          HeaderWrite( this, vec.size(), "VEC_PROBEIFC" );          /* WRITE */

          for ( auto & ifc : vec ) 
          {
            _write( os, ifc->pt );
            _write( os, ifc->props.sigtot );                              /* WRITE */
          }
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
        ASSERT ( vec.size() == h2.N , "Wrong vec size? " ) ; // << vec.size() << " != " << h2.N;
      }

      // QUADRATURE POINT LOOP
      for ( uint qp=0; qp<h2.N; qp++ )
      {
        RealTensor sigtot;         _read(is, sigtot);                      /* READ */
        if ( ! is_local ) continue;  // it is important to do the dummy runs above to keep sync
                                     //
        // Update the local structures
        vector<VPProps> & vec = by_elem.at(eid);
        VPProps &p = vec[qp];

        // TODO: currently we are only reading sigtot (file) => initial_stress (database)
        //       do we need to read other stuff ?
        p.initial_stress = sigtot;
      }
    }

    // The structure to feed
    ProbeByPnameByElemMap & probes_by_pname_by_elem = vp_ifc.probes_by_pname_by_elem;
    //

    /** Read Probes **/
    HeaderRead h3( this, "PROBES_BYBY" );              /* READ */

    for ( uint i=0; i<h3.N; i++ )
    {
      string pname; _read(is, pname);                  /* READ */

//      ASSERT( probes_by_pname_by_elem.count(pname) ) << "pname '" << pname << "'not found in gl_map.";
      auto & by_elem = probes_by_pname_by_elem[ pname ];  // Creates if non existing... No problem

      HeaderRead h4( this, "PROBES_BY_ELEM" );         /* READ */
      for ( uint j=0; j<h4.N; j++ )
      {
        uint eid; _read(is, eid);
        HeaderRead h5( this, "VEC_PROBEIFC" );         /* READ */
        bool is_local = ( mesh.processor_id() == mesh.elem_ref(eid).processor_id() );

        // Validation
//        if ( is_local ) 
//        {
//          if ( ! by_elem.count(eid) ) 
//          {
//            dlog(1) << "ELEMENTS IN PROCESSOR " << mesh.processor_id() << " (" << world.rank() << ") // pmame:" << pname << "  eid(" << eid << ") - active?" << mesh.elem_ref(eid).active();
//            for ( auto & [k,v] : by_elem ) dlog(1) << k;
//          }

//          ASSERT( by_elem.count(eid) , "Element not found in by_elem." );
//          auto & vec = by_elem.at( eid );

//          ASSERT( vec.size() == h5.N , "Wrong vec size? " ); // << vec.size() << " != " << h5.N;
//        }
        //
        if (is_local) 
        {
          if ( ! by_elem.count(eid) ) dlog(1) << " eid:" <<eid << "  N:" << h5.N;
          ASSERT( by_elem.count(eid) , "Element not found in by_elem." );
        }

        // For each point in the vector ...
        for ( uint pt_i=0; pt_i<h5.N; pt_i++ )
        {
          Point pt;          _read(is, pt);
          RealTensor sigtot; _read(is, sigtot);

          if ( ! is_local ) continue;

          if ( ! by_elem.count(eid) ) dlog(1) << "pt:" << pt << " eid:" <<eid;
          ASSERT( by_elem.count(eid) , "Element not found in by_elem." );

          auto & vec = by_elem.at( eid );
          ProbeIFC * ifc = vec[pt_i] ;
          
          ifc->props.initial_stress = sigtot;
        }
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
