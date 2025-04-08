
#include "material/ViscoPlasticMaterial.h"
#include "config/MaterialConfig.h"

#include "libmesh/system.h"
#include "libmesh/equation_systems.h"

#include "libmesh/string_to_enum.h"
#include "libmesh/enum_fe_family.h"
#include "libmesh/enum_order.h"

#include "libmesh/nonlinear_implicit_system.h"
#include "libmesh/transient_system.h"


/**
 *
 */
ViscoPlasticMaterial::ViscoPlasticMaterial( suint sid_, const MaterialConfig & config, System & sys_ ) :
  Material( sid_, config ), dof_indices_var(3),
  system( dynamic_cast<TransientNonlinearImplicitSystem&>(sys_) ),
  implicit(1),
  lame_mu( *(config.young) / 2 / ( 1 + *(config.poisson) ) ),
  lame_lambda( *(config.young) * *(config.poisson) / (1 + *(config.poisson)) / (1-2 * *(config.poisson)) ),
  bc_material(0)
{
  setup_variables();
}

ViscoPlasticMaterial::~ViscoPlasticMaterial()
{
  if ( bc_material ) delete(bc_material); 
  bc_material = 0;
}

/**
 *    Creates an identical Material, but prepared for integrating over boundaries.
 */
Material * ViscoPlasticMaterial::get_bc_material( Elem & elem, uint side, bool reinit )
{
  SCOPELOG(5);
  if ( ! bc_material ) 
  {
    bc_material = new ViscoPlasticMaterialBC( sid, config, system );
    bc_material->init_fem();
  }

  if ( reinit ) 
    bc_material->reinit( elem, side );

  return bc_material;
}

/**
 *
 */
ViscoPlasticMaterialBC::ViscoPlasticMaterialBC( suint sid_, const MaterialConfig & config_, System & sys_ ) :
  ViscoPlasticMaterial( sid_, config_, sys_ )
{
}

/**
 *
 */
void ViscoPlasticMaterial::setup_variables()
{
  // Add displacements variables to the current subdomain ID
  if (! config.fem_by_var.count( "U" ) ) flog << "Undefined var setup for variable 'U'. Please revise model material.FEM section.";
  auto & femspec = config.fem_by_var.at("U");
  Order order = Utility::string_to_enum<Order>( femspec.order ) ;
  FEFamily fe_family = Utility::string_to_enum<FEFamily>( femspec.family );

  dlog(1) << "Setting up variable 'U' for ViscoPlasticMaterial (sid=" << sid << ") ...";
  dlog(1) << "     Order:" << order;
  dlog(1) << "     FEFamily:" << fe_family;

  set<subdomain_id_type> sids = { sid };
  system.add_variable( "UX", order, fe_family, &sids );
  system.add_variable( "UY", order, fe_family, &sids );
  system.add_variable( "UZ", order, fe_family, &sids );

  implicit = femspec.implicit;
}

/**
 *  This can only be done after EquationSystems init.
 *  This is only done once for the material.
 */
void ViscoPlasticMaterial::init_fem()
{
  SCOPELOG(5);
  uint vid = system.variable_number( "UX" );

  // Setup shape functions
  DofMap & dof_map = system.get_dof_map();
  FEType fe_type = dof_map.variable_type(vid);

  fe = move( FEBase::build(3, fe_type) );

  // Setup gauss quadrature
  if ( ! is_bc() )
    qrule = QGauss( 3, fe_type.default_quadrature_order() );
  else 
    qrule = QGauss( 2, fe_type.default_quadrature_order() );

  fe->attach_quadrature_rule (&qrule);

  // Enable calculations calculations
  fe->get_JxW();
  fe->get_phi();
  fe->get_dphi();

  if ( is_bc() ) 
    fe->get_normals(); 

  // Jacobian
  for ( uint i=0; i<3; i++ )
  {
    vector<DenseSubMatrix<Number>> kk;
    for ( uint j=0; j<3; j++ ) kk.emplace_back( Ke );
    Ke_var.push_back(kk);
  }

  // RHS Vector
  for ( uint i=0; i<3; i++ )
    Re_var.emplace_back(Re);
}

/**
 *  Init the DoF map and the element matrices.
 *  This function must be called before (or at the beginning of)
 *  the jacobian and residual.
 *
 *  _side_ is an optional parameter, used only when the material
 *  is being initialized for a BC.
 */
void ViscoPlasticMaterial::reinit( const Elem & elem, uint side )
{
  SCOPELOG(5);

  if ( is_bc() ) 
    fe->reinit( &elem, side );
  else
    fe->reinit( &elem );
  
  const DofMap & dof_map = system.get_dof_map();

  dof_map.dof_indices (&elem, dof_indices);
  dof_map.dof_indices (&elem, dof_indices_var[0], 0);
  dof_map.dof_indices (&elem, dof_indices_var[1], 1);
  dof_map.dof_indices (&elem, dof_indices_var[2], 2);

  uint n_dofs = dof_indices.size();
  uint n_dofsv = dof_indices_var[0].size();
  Ke.resize (n_dofs, n_dofs);

  for (uint vi=0; vi<3; vi++)
  for (uint vj=0; vj<3; vj++)
    Ke_var[vi][vj].reposition ( vi*n_dofsv, vj*n_dofsv, n_dofsv, n_dofsv );

  Re.resize (n_dofs);
  for (uint vi=0; vi<3; vi++)
    Re_var[vi].reposition ( vi*n_dofsv, n_dofsv );
}

/**
 *     Builds the jacobian of the element and assembles in the global _jacobian_.
 */
void ViscoPlasticMaterial::jacobian (const NumericVector<Number> & soln, SparseMatrix<Number> & jacobian)
{
  SCOPELOG(5);
  const vector<Real> & JxW = fe->get_JxW();
  const vector<vector<Real>> & phi = fe->get_phi();
  const vector<vector<RealGradient>> & dphi = fe->get_dphi();
  const DofMap & dof_map = system.get_dof_map();
  uint n_dofsv = dof_indices_var[0].size();

  // ****
  // Effective mechanics
  //       ( \phi_i,j , C_ijkl u_k,l )   ok
  for (uint qp=0; qp<qrule.n_points(); qp++   )
  for (uint B=0; B<n_dofsv;  B++)
  for (uint M=0; M<n_dofsv;  M++)
  for (uint i=0; i<3; i++) 
  for (uint j=0; j<3; j++)
  for (uint k=0; k<3; k++) 
  for (uint l=0; l<3; l++) 
    Ke_var[i][k](B,M) += JxW[qp] * C_ijkl(i,j,k,l) * dphi[M][qp](l) * dphi[B][qp](j);
//    Ke_var[i][k](B,M) += implicit * JxW[qp] * C_ijkl(i,j,k,l) * dphi[M][qp](l) * dphi[B][qp](j);

//  dlog(1) << endl << Ke;
  // Add the the global matrix
  dof_map.constrain_element_matrix (Ke, dof_indices);
  jacobian.add_matrix (Ke, dof_indices);
}



/**
 *     Builds the RHS of the element and assembles in the global _residual_.
 */
void ViscoPlasticMaterial::residual (const NumericVector<Number> & soln, NumericVector<Number> & residual)
{
  SCOPELOG(5);
  const vector<Real> & JxW = fe->get_JxW();
  const vector<vector<Real>> & phi = fe->get_phi();
  const vector<vector<RealGradient>> & dphi = fe->get_dphi();
  const DofMap & dof_map = system.get_dof_map();
  uint n_dofsv = dof_indices_var[0].size();

  // Computes grad_u 
  vector<TensorValue<Number>> OLD_GRAD_U(qrule.n_points());
  vector<TensorValue<Number>> GRAD_U(qrule.n_points());
  for (uint qp=0;    qp<qrule.n_points(); qp++   )
  for (uint i=0; i<3; i++)
  for (uint j=0; j<3; j++)
  for (uint M=0;  M<n_dofsv;  M++)
  {
    OLD_GRAD_U[qp](i, j) += dphi[M][qp](j) * system.old_solution(dof_indices_var[i][M]);
    GRAD_U[qp](i, j) += dphi[M][qp](j) * soln(dof_indices_var[i][M]);
  }

  // ( \phi_j , Cijkl U_k,l ) 
  for (uint qp=0; qp<qrule.n_points(); qp++   )
  for (uint B=0;  B<n_dofsv;  B++)
  for (uint i=0; i<3; i++) 
  for (uint j=0; j<3; j++) 
  for (uint k=0; k<3; k++) 
  for (uint l=0; l<3; l++) 
  {
    Re_var[i](B) += JxW[qp] *  dphi[B][qp](j) * C_ijkl(i,j,k,l) * GRAD_U[qp](k,l) ;
//    Re_var[i](B) += implicit     * JxW[qp] *  dphi[B][qp](j) * C_ijkl(i,j,k,l) * GRAD_U[qp](k,l) ;
//    Re_var[i](B) -= (1-implicit) * JxW[qp] *  dphi[B][qp](j) * C_ijkl(i,j,k,l) * OLD_GRAD_U[qp](k,l) ;
  }

  dof_map.constrain_element_vector ( Re, dof_indices );
  residual.add_vector ( Re, dof_indices );
}


/**
 *
 *
 *
 *  BOUNDARY CONSTRAINTS
 *
 *
 *
 */


/**
 *     Builds the jacobian of the element and assembles in the global _jacobian_.
 */
void ViscoPlasticMaterialBC::jacobian (const NumericVector<Number> & soln, SparseMatrix<Number> & jacobian)
{
  SCOPELOG(5);
  const vector<Real> & JxW = fe->get_JxW();
  const vector<vector<Real>> & phi = fe->get_phi();
  const vector<vector<RealGradient>> & dphi = fe->get_dphi();
  const DofMap & dof_map = system.get_dof_map();
  uint n_dofsv = dof_indices_var[0].size();

  // Nothing to do here.
}

/**
 *     Builds the RHS of the element and assembles in the global _residual_.
 */
void ViscoPlasticMaterialBC::residual (const NumericVector<Number> & soln, NumericVector<Number> & residual)
{
  SCOPELOG(5);
  const vector<Real> & JxW = fe->get_JxW();
  const vector<vector<Real>> & phi = fe->get_phi();
  const vector<vector<RealGradient>> & dphi = fe->get_dphi();
  const DofMap & dof_map = system.get_dof_map();
  const vector<Point> & normals = fe->get_normals(); 
  uint n_dofsv = dof_indices_var[0].size();

  if ( sigtot )
  for (uint qp=0; qp<qrule.n_points(); qp++) 
  for (uint B=0; B<n_dofsv; B++)
  for (uint i=0; i<3; i++) 
  for (uint j=0; j<3; j++) 
    Re_var[i](B) += JxW[qp] * (*sigtot)(i,j) * normals[qp](j) * phi[B][qp];
//    Re_var[i](B) += (1-implicit) * JxW_face[qp] * (*sigtot)(i,j) * normals[qp](j) * phi_u_face[B][qp];

  dof_map.constrain_element_vector (Re, dof_indices);
  residual.add_vector (Re, dof_indices);
}


