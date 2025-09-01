
#include "solver/viscoplastic/EGDataStruct.h"

#include "libmesh/dof_map.h"

namespace solver {
namespace viscoplastic {

EGFEM::EGFEM( System & sys ) 
{
  if ( ! sys.has_variable( "UegX" ) )    // We are not EG?
    flog << "Only EG systems should be here.";

  uint vid = sys.variable_number( "UegX" );

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
void EGFEM::attach_qrule( QBase * qrule ) {
  fe->attach_quadrature_rule( qrule ); 
}

/** **/
void EGFEM::set_dofs( System & sys, const Elem * elem )
{
  const DofMap & dof_map = sys.get_dof_map();
  vector<dof_id_type> di;

  dofi_cg.clear();
  dofi_eg.clear();

  // 0..5 = "UX", "UY", "UZ", "UegX", "UegY", "UegZ"

  // CG
  for ( uint vi=0; vi<3; vi++ )
  {
    dof_map.dof_indices ( elem, di, vi );
    dofi_cg.insert( dofi_cg.end(), di.begin(), di.end() );
  }

  // EG
  for ( uint vi=3; vi<6; vi++ )
  {
    dof_map.dof_indices ( elem, di, vi );
    dofi_eg.insert( dofi_eg.end(), di.begin(), di.end() );
  }

}

/** **/
ostream& operator<<(ostream& os, const EGFace & m)
{
  os << endl;
  os << "       EGFace" << endl;
  os << "            E:S : " << m.eid << ":" << m.side << endl; 
  os << "            Nqp : " << m.Pq.size() << endl;

  return os;
}

/** **/
ostream& operator<<(ostream& os, const EGFacePair & m)
{
  os << endl;
  os << "       EGFacePair" << endl;
  os << "            E:S+ : " << m.eid_p << ":" << m.side_p << endl; 
  os << "            E-   : " << m.eid_n << endl;
  os << "            Nqp  : " << m.Pq_p.size() << " , " << m.Pq_n.size() << endl;
  return os;
}

}} // ns
