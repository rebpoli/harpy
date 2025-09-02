
#include "solver/viscoplastic/EGDataStruct.h"
#include "util/OutputOperators.h"

#include "libmesh/dof_map.h"

namespace solver {
namespace viscoplastic {

EGFEM::EGFEM( System & sys ) 
{
  if ( ! sys.has_variable( "UegX" ) )    // We are not EG?
    flog << "Only EG systems should be here.";

  // EG
  {
    uint vid = sys.variable_number( "UegX" );
    // Setup shape functions
    DofMap & dof_map = sys.get_dof_map();
    FEType fe_type = dof_map.variable_type(vid);
    fe_eg = move( FEBase::build(3, fe_type) );
    // Enable calculations calculations
    fe_eg->get_JxW();
    fe_eg->get_phi();
    fe_eg->get_dphi();
    fe_eg->get_xyz();
    fe_eg->get_normals();
  }

  // CG
  {
    uint vid = sys.variable_number( "UX" );
    // Setup shape functions
    DofMap & dof_map = sys.get_dof_map();
    FEType fe_type = dof_map.variable_type(vid);
    fe_cg = move( FEBase::build(3, fe_type) );
    // Enable calculations calculations
    fe_cg->get_JxW();
    fe_cg->get_phi();
    fe_cg->get_dphi();
    fe_cg->get_xyz();
    fe_cg->get_normals();
  }
}

/** **/
void EGFEM::attach_qrule( QBase * qrule ) {
  fe_cg->attach_quadrature_rule( qrule ); 
  fe_eg->attach_quadrature_rule( qrule ); 
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
  using util::operator<<;
  os << endl;
  os << "       EGFacePair" << endl;
  os << "            E:S+ : " << m.eid_p << ":" << m.side_p << endl; 
  os << "            E-   : " << m.eid_n << endl;
  os << "            Nqp  : " << m.Pq_p.size() << " , " << m.Pq_n.size() << endl;
  os << "            Pq_p  : " << m.Pq_p << endl;
  os << "            Pq_n  : " << m.Pq_n << endl;
  return os;
}

}} // ns
