#pragma once

#include "harpy/Global.h"
#include "solver/viscoplastic/VPProps.h"

#include "libmesh/system.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/fe.h"

#include <optional>

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

  inline VPProps * get_P( uint qp ) {
    ASSERT( qp < Pq.size() , "Out of bounds request in EGFacePair::get_Pp. QP=" << qp << " size:" << Pq.size() );
    return & ( Pq.at(qp) ); 
  }
};

/**
 *   The description of a boundary, for a variable,
 *   for a value
 */
struct EGDirichlet 
{
  string vname;         // the CG name (UX,UY,UZ)
  uint vid_cg, vid_eg;
  double val;
  vector<EGFace> egface_vec;

  inline vector<optional<double>> u_hat()
  {
    vector<optional<double>> ret(3);
    if ( vname == "UX" ) ret[0] = val;
    if ( vname == "UY" ) ret[1] = val;
    if ( vname == "UZ" ) ret[2] = val;
    return ret;
  }

  inline string vname_eg()
  {
    if ( vname == "UX" ) return "UegX";
    if ( vname == "UY" ) return "UegY";
    if ( vname == "UZ" ) return "UegZ";
    flog << "Unknonw variable for Dirichlet setting: " << vname << ".";
    return vname;
  }
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
  unique_ptr<FEBase> fe_cg, fe_eg;               /// The finite element object to hold shape funtions, jxw, etc
  vector<dof_id_type> dofi_eg, dofi_cg;

  void attach_qrule( QBase * qrule );
  void set_dofs( System & sys, const Elem * elem );
};

/** **/
ostream& operator<<(ostream& os, const EGFace & m);
ostream& operator<<(ostream& os, const EGFacePair & m);

}} // ns
