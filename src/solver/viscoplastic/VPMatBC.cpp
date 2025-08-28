#include "solver/viscoplastic/VPMatBC.h"
#include "config/MaterialConfig.h"
#include "solver/viscoplastic/VPSolver.h"

#include "libmesh/fe.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/dof_map.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/dense_submatrix.h"
#include "libmesh/dense_subvector.h"

namespace solver {
namespace viscoplastic {

/**
 *
 */
ViscoPlasticMaterialBC::ViscoPlasticMaterialBC( suint sid_, const MaterialConfig * config_,
                                                TransientNonlinearImplicitSystem & sys_,
                                                ViscoplasticSolver & vpsolver_ ) :
  config( config_ ),
  QP(0),
  vpsolver(vpsolver_), 
  system( sys_ ), 
  qrule(3),
  elem(0),
  n_dofs(0), n_dofsv(0), n_uvars( vpsolver_.is_eg() ? 6 : 3 ), n_dofs_eg(0),
  name(config->name),
  sid( sid_ )
{
}

/**
 *  This can only be done after EquationSystems init.
 *  This is only done once for the material.
 */
void ViscoPlasticMaterialBC::init_fem()
{
  SCOPELOG(1);
  uint vid = system.variable_number( "UX" );

  // Setup shape functions
  DofMap & dof_map = system.get_dof_map();
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

/**
 *
 */
void ViscoPlasticMaterialBC::reinit( const NumericVector<Number> & soln, const Elem & elem_, uint side )
{
  SCOPELOG(5);


  QP = 0;
  elem = &elem_;

  fe->reinit( elem, side );
  
  const DofMap & dof_map = system.get_dof_map();
  dof_map.dof_indices (elem, dof_indices);

  vector<dof_id_type> dof_indices_var;
  dof_map.dof_indices (elem, dof_indices_var, 0);

  n_dofs = dof_indices.size();
  n_dofsv = dof_indices_var.size();

  if ( vpsolver.is_eg() ) // EG
  {
    vector<dof_id_type> dof_indices_eg;
    dof_map.dof_indices (elem, dof_indices_eg, 3);
    n_dofs_eg = dof_indices_eg.size();
    dlog(1) << "n_dofs_eg: " << n_dofs_eg;
  }

  Re.resize (n_dofs);

  // Update the current element in the interfacse
//  uint eid = elem->id();
//  uint nqp = qrule.size();

  next_qp(0);

  // Initialize the flattened containers
  _init_autodiff();

  // Feed Uib with current solution
  for ( uint i=0; i<n_uvars; i++ )
  for ( uint B=0; B<Ndof(i); B++ )
    Uib(i,B) = soln( dof_indices[_idx(i,B)] );

}

/**
 *
 */
void ViscoPlasticMaterialBC::residual_and_jacobian_qp ()
{
  SCOPELOG(5);
  const vector<Real> & JxW = fe->get_JxW();
  const vector<vector<Real>> & phi = fe->get_phi();
  const vector<Point> & normals = fe->get_normals(); 

  ad_Fib.setZero(); 

  // Note: as the sigtot is constant, we dont need the jacobian
  //     : Keeping the variable as a AUTODIFF for parallelism
  //     : with the parent object. 
  if ( sigtot )
  for (uint B=0; B<n_dofsv; B++)
  for (uint i=0; i<n_uvars; i++) 
  for (uint j=0; j<n_uvars; j++) 
    Fib(i,B) -= JxW[QP] * (*sigtot)(i,j) * normals[QP](j) * phi[B][QP];

  // Map from the AD variable to libmesh datastructures
  for (uint i=0; i<n_uvars; i++) 
  for (uint B=0;  B<n_dofsv;  B++)
    Re( i*n_dofsv + B ) += val( Fib(i,B) );

}

/**
 *     Builds the RHS of the element and assembles in the global _residual_.
 */
void ViscoPlasticMaterialBC::residual_and_jacobian ( Elem & elem_, uint side, 
                                                     const NumericVector<Number> & soln, 
                                                     SparseMatrix<Number> * jacobian ,
                                                     NumericVector<Number> * residual )
{
  UNUSED(jacobian);
  // Reinit the object
  reinit( soln, elem_, side );

  // Build the element jacobian and residual for each quadrature point _qp_.
  do { residual_and_jacobian_qp(); } while ( next_qp() );

  // Add to the global residual vector
  if ( residual )
  {
    const DofMap & dof_map = system.get_dof_map();
    dof_map.constrain_element_vector (Re, dof_indices);
    residual->add_vector (Re, dof_indices);
  }
}



}} // ns
