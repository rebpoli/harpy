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
struct EGFace { 
  EGFace( uint e, uint s ) : eid(e), side(s) {}
  uint eid, side;
  vector< VPProps > Pq;         /// VPProps for each quadrature point
};

/**
 *
 */
struct EGFacePair { 
  EGFacePair( uint ep, uint sp, uint en ) : eid_p(ep), side_p(sp), eid_n(en) {}
  uint eid_p, side_p;
  uint eid_n;
  vector< VPProps > Pq_p;       /// VPProps for each quadrature point
  vector< VPProps > Pq_n;       /// VPProps for each quadrature point
};

/**
 *
 */
struct EGFEM {
  EGFEM( System & sys );
  const Elem * elem;                               /// The element that the material has been reinit'ed to
  vector<dof_id_type> dof_indices;
  unique_ptr<FEBase> fe;                           /// The finite element object to hold shape funtions, jxw, etc

  void attach_qrule( QBase * qrule );
};

/** **/
ostream& operator<<(ostream& os, const EGFace & m);
ostream& operator<<(ostream& os, const EGFacePair & m);

}} // ns
