#pragma once

#include "harpy/Global.h"
#include "solver/viscoplastic/VPProps.h"

#include "libmesh/system.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/fe.h"

namespace libMesh { class Elem; } 

namespace solver {
namespace viscoplastic {

using namespace libMesh;

/**
 *
 */
struct EGFace 
{ 
  EGFace( uint e, uint s ) : eid(e), side(s) {}
  uint eid, side;
  vector< VPProps > Pq;         /// VPProps for each quadrature point
};

/**
 *
 */
struct EGFacePair 
{ 
  EGFacePair( uint ep, uint sp, uint en ) : eid_p(ep), side_p(sp), eid_n(en) {}
  uint eid_p, side_p;
  uint eid_n;
  vector< VPProps > Pq_p;       /// VPProps for each quadrature point
  vector< VPProps > Pq_n;       /// VPProps for each quadrature point

  ///
  inline VPProps * get_Pp( uint qp ) {
    ASSERT( qp < Pq_p.size() , "Out of bounds request in EGFacePair::get_Pp. QP=" << qp << " size:" << Pq_p.size() );
    return & ( Pq_p.at(qp) ); 
  }
  ///
  inline VPProps * get_Pn( uint qp ) {
    ASSERT( qp < Pq_n.size() , "Out of bounds request in EGFacePair::get_Pn. QP=" << qp << "size:" << Pq_n.size() );
    return & ( Pq_n.at(qp) ); 
  }
};

/**
 *
 */
struct EGFEM 
{
  EGFEM( System & sys );
  unique_ptr<FEBase> fe;                           /// The finite element object to hold shape funtions, jxw, etc
  vector<dof_id_type> dofi_eg, dofi_cg;

  void attach_qrule( QBase * qrule );
  void set_dofs( System & sys, const Elem * elem );
};

/** **/
ostream& operator<<(ostream& os, const EGFace & m);
ostream& operator<<(ostream& os, const EGFacePair & m);

}} // ns
