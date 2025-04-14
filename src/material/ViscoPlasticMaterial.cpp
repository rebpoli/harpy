
#include "material/ViscoPlasticMaterial.h"
#include "config/MaterialConfig.h"

#include "libmesh/system.h"
#include "libmesh/equation_systems.h"

#include "libmesh/string_to_enum.h"
#include "libmesh/enum_fe_family.h"
#include "libmesh/enum_order.h"

#include "libmesh/nonlinear_implicit_system.h"
#include "libmesh/transient_system.h"

#include "util/OutputOperators.h"


/**
 *
 */
ViscoPlasticMaterial::ViscoPlasticMaterial( suint sid_,
                                            const MaterialConfig & config, 
                                            TransientNonlinearImplicitSystem & sys_ ) :
  Material( sid_, config ), dof_indices_var(3),
  system( sys_ ), implicit(1), bc_material(0)
{
  // Lists the necessary properties to fetch from the config during init_coupler
  required_material_properties.assign({
      "alpha_d", "beta_e", "young", "poisson", "lame_mu", "lame_lambda", 
  });

  // Setup libmesh system variables
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
Material * ViscoPlasticMaterial::get_bc_material()
{
  SCOPELOG(5);
  if ( ! bc_material ) 
  {
    bc_material = new ViscoPlasticMaterialBC( sid, config, system );
    bc_material->init_fem();
  }

  return bc_material;
}

/**
 *
 */
ViscoPlasticMaterialBC::ViscoPlasticMaterialBC( suint sid_, const MaterialConfig & config_,
                                                TransientNonlinearImplicitSystem & sys_ ) :
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
  fe->get_xyz();

  if ( is_bc() ) 
    fe->get_normals(); 

  // Jacobian
  for ( uint i=0; i<3; i++ )
  {
    vector<DenseSubMatrix<Number>> kk;
    for ( uint j=0; j<3; j++ ) kk.emplace_back( Ke );
    Ke_var.push_back(kk);
  }

//  // RHS Vector
//  for ( uint i=0; i<3; i++ )
//    Re_var.emplace_back(Re);

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

  elem_coupler = 0;

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
//  for (uint vi=0; vi<3; vi++)
//    Re_var[vi].reposition ( vi*n_dofsv, n_dofsv );
}

/**
 *     Calculates the output information at a single point in this material.
 *     Pushes the results into the trg_coupler
 *     _elem_ is the element in this mesh where the point lies
 */
void ViscoPlasticMaterial::feed_coupler( ElemCoupler & trg_ec, const Point & trg_pt,
                                         const Elem * elem, const NumericVector<Number> & soln )
{
  SCOPELOG(5);

  // Interpolators are the ones from this material.
  // However, they must be initialized with the xyz of the gauss points of the target
  const vector<vector<Real>> & phi = fe->get_phi();
  const vector<vector<RealGradient>> & dphi = fe->get_dphi();
  const DofMap & dof_map = system.get_dof_map();
  dof_map.dof_indices (elem, dof_indices);
  dof_map.dof_indices (elem, dof_indices_var[0], 0);
  uint n_dofsv = dof_indices_var[0].size();

  std::vector<Point> pts_from = { trg_pt };
  std::vector<Point> pts_to;
  FEInterface::inverse_map (3, fe->get_fe_type(), elem, pts_from, pts_to, 1E-10);
  fe->reinit( elem, & pts_to );

  // Prepare the Uib vector for the automatic differentiation
  Uib.clear();
  for ( uint i=0; i<3; i++ )
  {
    vector<Number> row;
    for ( uint B=0; B<n_dofsv; B++ )
      row.push_back( soln( dof_indices[i*n_dofsv + B] ) );
    Uib.push_back(row);
  }

  RealVectorValue U = 0;
  for ( uint B=0; B<n_dofsv; B++ )
  for ( uint i=0; i<3; i++ )
    U(i) += phi[B][0] * Uib[i][B];

  RealTensor GRAD;
  for ( uint B=0; B<n_dofsv; B++ )
  for ( uint i=0; i<3; i++ )
  for ( uint j=0; j<3; j++ )
    GRAD(i,j) += dphi[B][0](j) * Uib[i][B];

  // Feed the coupler in this point
  vector<RealVectorValue> & trg_Uqi  = trg_ec.vector_params["U"];
  trg_Uqi.push_back( U );

  vector<RealTensor> & trg_GRAD_Uqij = trg_ec.tensor_params["GRAD_U"];
  trg_GRAD_Uqij.push_back( GRAD );
}

/**
 *  Init the DoF map and the element matrices.
 *
 *  Then, initializes the solution depending on soln.
 */
void ViscoPlasticMaterial::reinit( const NumericVector<Number> & soln, const Coupler & coupler, const Elem & elem, uint side )
{
  SCOPELOG(5);

  /*
   * Initializes the FEM structures, not depending on the solution
   */
  reinit( elem, side );

  /*
   * Now initializes the structures depending on the solution vector
   */
  uint n_dofsv = dof_indices_var[0].size();

  // Prepare the Uib vector for the automatic differentiation
  Uib.clear();
  for ( uint i=0; i<3; i++ )
  {
    vector<Number> row;
    for ( uint B=0; B<n_dofsv; B++ )
      row.push_back( soln( dof_indices[i*n_dofsv + B] ) );
    Uib.push_back(row);
  }

  // Prepare Fib vector for the automatic differentiation
  Fib.clear();
  for ( uint i=0; i<3; i++ )
  {
    vector<Number> row;
    for ( uint B=0; B<n_dofsv; B++ ) row.push_back( 0 );
    Fib.push_back(row);
  }

  uint eid = elem.id();
  if ( ! coupler.count( eid ) ) flog << "Coupler has not been initialized for element '" << eid << "'. Something is terribly wrong.";
  elem_coupler = & ( coupler.at(eid) );

  // Fetch the needed parameters from the coupler
  get_from_element_coupler(  "T",             temperature    ); 
  get_from_element_coupler(  "alpha_d",       alpha_d        ); 
  get_from_element_coupler(  "lame_mu",       lame_mu        ); 
  get_from_element_coupler(  "lame_lambda",   lame_lambda    ); 
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
    Ke_var[i][k](B,M) += JxW[qp] * C_ijkl(qp, i,j,k,l) * dphi[M][qp](l) * dphi[B][qp](j);
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
    GRAD_U[qp](i, j) += dphi[M][qp](j) * Uib[i][M];

  // ( \phi_j , Cijkl U_k,l ) 
  for (uint qp=0; qp<qrule.n_points(); qp++   )
  for (uint B=0;  B<n_dofsv;  B++)
  for (uint i=0; i<3; i++) 
  for (uint j=0; j<3; j++) 
  for (uint k=0; k<3; k++) 
  for (uint l=0; l<3; l++) 
    Fib[i][B] += JxW[qp] *  dphi[B][qp](j) * C_ijkl(qp, i,j,k,l) * GRAD_U[qp](k,l) ;

  // Temperature
  //       alpha_d = beta_d * K_bulk
  //       - ( \phi_i,i , \alpha_d T ) ==> term in the RHS.
  for (uint qp=0; qp<qrule.n_points(); qp++   )
  for (uint B=0;  B<n_dofsv;  B++)
  for (uint i=0;  i<3;  i++)
    Fib[i][B] -= JxW[qp] *  dphi[B][qp](i) * alpha_d[qp]  * temperature[qp];

  /*
  *    Build the Re vector and assemble in the global residual vector
  */
  for (uint i=0; i<3; i++) 
  for (uint B=0;  B<n_dofsv;  B++)
    Re( i*n_dofsv + B ) = Fib[i][B];

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
    Fib[i][B] -= JxW[qp] * (*sigtot)(i,j) * normals[qp](j) * phi[B][qp];
//    Re_var[i](B) += (1-implicit) * JxW_face[qp] * (*sigtot)(i,j) * normals[qp](j) * phi_u_face[B][qp];

  // Build the Re vector
  for (uint i=0; i<3; i++) 
  for (uint B=0;  B<n_dofsv;  B++)
    Re( i*n_dofsv + B ) = Fib[i][B];

  dof_map.constrain_element_vector (Re, dof_indices);
  residual.add_vector (Re, dof_indices);
}


