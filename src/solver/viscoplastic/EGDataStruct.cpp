
#include "solver/viscoplastic/EGDataStruct.h"

#include "libmesh/dof_map.h"

namespace solver {
namespace viscoplastic {

EGFEM::EGFEM( System & sys ) : 
  qrule(2),
  elem(0) 
{
  if ( ! sys.has_variable( "UegX" ) )    // We are not EG?
    flog << "Only EG systems should be here.";

  uint vid = sys.variable_number( "UegX" );

  // Setup shape functions
  DofMap & dof_map = sys.get_dof_map();
  FEType fe_type = dof_map.variable_type(vid);

  fe = move( FEBase::build(3, fe_type) );

  // Setup gauss quadrature
  qrule = QGauss( 2, fe_type.default_quadrature_order() );

  fe->attach_quadrature_rule (&qrule);

  // Enable calculations calculations
  fe->get_JxW();
  fe->get_phi();
  fe->get_dphi();
  fe->get_xyz();
  fe->get_normals();
}


}} // ns
