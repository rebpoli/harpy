
#include "solver/viscoplastic/DGDataStruct.h"
#include "util/OutputOperators.h"

#include "libmesh/dof_map.h"

namespace solver {
namespace viscoplastic {

DGFEM::DGFEM( System & sys ) 
{
  uint vid = sys.variable_number( "UX" );
  // Setup shape functions
  DofMap & dof_map = sys.get_dof_map();
  FEType fe_type = dof_map.variable_type(vid);
  fe = move( FEBase::build(3, fe_type) );
  // Enable calculations calculations
  fe->get_JxW();
  fe->get_phi();
  fe->get_dphi();
  fe->get_xyz();
  fe->get_normals();
}

/** **/
void DGFEM::attach_qrule( QBase * qrule ) {
  fe->attach_quadrature_rule( qrule ); 
}

/** **/
void DGFEM::set_dofs( System & sys, const Elem * elem )
{
  const DofMap & dof_map = sys.get_dof_map();
  vector<dof_id_type> di;

  dofi.clear();

  for ( uint vi=0; vi<3; vi++ )
  {
    dof_map.dof_indices ( elem, di, vi );
    dofi_cg.insert( dofi_cg.end(), di.begin(), di.end() );
  }
}

/** **/
ostream& operator<<(ostream& os, const DGFace & m)
{
  os << m.eid << "/" << m.side << "(nqp:" << m.Pq.size() << ") "; 
  return os;
}

/** **/
ostream& operator<<(ostream& os, const DGFacePair & m)
{
  using util::operator<<;
  os << endl;
  os << "       DGFacePair" << endl;
  os << "            E:S+ : " << m.eid_p << ":" << m.side_p << endl; 
  os << "            E-   : " << m.eid_n << endl;
  os << "            Nqp  : " << m.Pq_p.size() << " , " << m.Pq_n.size() << endl;
  os << "            Pq_p  : " << m.Pq_p << endl;
  os << "            Pq_n  : " << m.Pq_n << endl;
  return os;
}

}} // ns
