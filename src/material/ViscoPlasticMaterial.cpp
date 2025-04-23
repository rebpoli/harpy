
#include "material/ViscoPlasticMaterial.h"
#include "config/MaterialConfig.h"
#include "solver/ViscoplasticSolver.h"

#include "libmesh/system.h"
#include "libmesh/equation_systems.h"

#include "libmesh/string_to_enum.h"
#include "libmesh/enum_fe_family.h"
#include "libmesh/enum_order.h"

#include "libmesh/nonlinear_implicit_system.h"
#include "libmesh/transient_system.h"

#include "util/OutputOperators.h"
#include "harpy/Timestep.h"


/**
 *    BC_OBJECT is true when this is being called from the children constructor
 */
ViscoPlasticMaterial::ViscoPlasticMaterial( suint sid_,
                                            const MaterialConfig & config, 
                                            TransientNonlinearImplicitSystem & sys_,
                                            ViscoplasticSolver & vpsolver_,
                                            bool called_from_bc_constructor ) :
  Material( sid_, config ),
  dof_indices_var(3),
  P( 0 ), T(0), bc_material(0),
  vpsolver(vpsolver_), system( sys_ ), 
  stress_system( vpsolver.stress_system ),
  stress_postproc( *this, stress_system )
{
  // Setup cariables only if on the parent
  if (! called_from_bc_constructor ) 
  {
    setup_variables(); /// Only the parent material need to setup the variables

    // Setup the BC material only on the parent
    bc_material = new ViscoPlasticMaterialBC( sid, config, system, vpsolver );
  }
}

/**
 *
 */
ViscoPlasticMaterial::~ViscoPlasticMaterial()
{
  if ( bc_material ) delete(bc_material); 
  bc_material = 0;
}

/**
 *    This function is called every reinit() of the material.
 *    It updates the internal property pointer (*OP) 
 *
 *    The inteface must have been initialized.
 */
void ViscoPlasticMaterial::init_properties()
{
  SCOPELOG(5);
  // Calc xyz (fe must have already been initialized)
  const std::vector<Point> & xyz = fe->get_xyz();

  // We have it in cache? Return;
  if ( vp_ifc.valid ) return; 

  uint qp = 0;
  for ( auto & pt : xyz )
  {
    auto & prop = vp_ifc.get( qp );

    prop.alpha_d          = config.get_property( "alpha_d",         pt,     "porothermoelastic" );
    prop.beta_e           = config.get_property( "beta_e",          pt,     "porothermoelastic" );
    prop.lame_mu          = config.get_property( "lame_mu",         pt,     "porothermoelastic" );
    prop.lame_lambda      = config.get_property( "lame_lambda",     pt,     "porothermoelastic" );
    prop.creep_carter_a   = config.get_property( "a",               pt,     "creep_carter" );
    prop.creep_carter_q   = config.get_property( "q",               pt,     "creep_carter" );
    prop.creep_carter_n   = config.get_property( "n",               pt,     "creep_carter" );
    qp++;
  }
  vp_ifc.valid = 1;
}


/**
 *    Creates an identical Material, but prepared for integrating over boundaries.
 */
ViscoPlasticMaterialBC * ViscoPlasticMaterial::get_bc_material()
{
  SCOPELOG(5);
  if ( ! bc_material ) flog << "BC material not set? It should have been set in the constructor of the parent class, for all materials.";

  return bc_material;
}

/**
 *
 */
ViscoPlasticMaterialBC::ViscoPlasticMaterialBC( suint sid_, const MaterialConfig & config_,
                                                TransientNonlinearImplicitSystem & sys_,
                                                ViscoplasticSolver & vpsolver_ ) :
  ViscoPlasticMaterial( sid_, config_, sys_, vpsolver_, true )
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

  stress_postproc.init_fem();

  if ( bc_material ) bc_material->init_fem();
}

/**
 *  Init the DoF map and the element matrices.
 *  This function must be called before (or at the beginning of)
 *  the jacobian and residual.
 *
 *  _side_ is an optional parameter, used only when the material
 *  is being initialized for a BC.
 */
void ViscoPlasticMaterial::reinit( const Elem & elem_, uint side )
{
  SCOPELOG(5);

  QP = 0;
  elem = &elem_;

  if ( is_bc() ) 
    fe->reinit( elem, side );
  else
    fe->reinit( elem );
  
  const DofMap & dof_map = system.get_dof_map();

  dof_map.dof_indices (elem, dof_indices);
  dof_map.dof_indices (elem, dof_indices_var[0], 0);
  dof_map.dof_indices (elem, dof_indices_var[1], 1);
  dof_map.dof_indices (elem, dof_indices_var[2], 2);

  uint n_dofs = dof_indices.size();
  uint n_dofsv = dof_indices_var[0].size();
  Ke.resize (n_dofs, n_dofs);

  for (uint vi=0; vi<3; vi++)
  for (uint vj=0; vj<3; vj++)
    Ke_var[vi][vj].reposition ( vi*n_dofsv, vj*n_dofsv, n_dofsv, n_dofsv );

  Re.resize (n_dofs);

  // Update the current element in the interfacse
  uint eid = elem->id();
  uint nqp = qrule.size();

  // Initialize interfacse
  vp_ifc.reinit( eid, nqp );
  th_ifc.reinit( eid, nqp );
  next_qp(0);

  /// Init properties from confguration file if needed. Updates material_properties_by_qp
  init_properties();

}

/**
 *  Init the DoF map and the element matrices.
 *
 *  Then, initializes the solution depending on soln.
 */
void ViscoPlasticMaterial::reinit( const NumericVector<Number> & soln, const Elem & elem_, uint side )
{
  SCOPELOG(5);

  /*
   * Initializes the FEM structures, not depending on the solution
   */
  reinit( elem_, side );

  /*
   * Now initializes the structures depending on the solution vector
   */
  uint n_dofsv = dof_indices_var[0].size();

  // Initialize the flattened containers
  _init_autodiff( n_dofsv );

  // Feed Uib with current solution
  for ( uint i=0; i<3; i++ )
  for ( uint B=0; B<n_dofsv; B++ )
    Uib(i,B) = soln( dof_indices[i*n_dofsv + B] );

}

/**
 *     Updates the interface based on the current Uib solution
 */
void ViscoPlasticMaterial::update_ifc_qp()
{
//  SCOPELOG(10);
//  uint n_dofsv = dof_indices_var[0].size();
//  const vector<vector<RealGradient>> & dphi = fe->get_dphi();

//  // Computes grad_u  ==> move to reinit, and to the interface (?)
//  P->grad_u = 0;
//  for (uint i=0; i<3; i++)
//  for (uint j=0; j<3; j++)
//  for (uint M=0;  M<n_dofsv;  M++)
//    P->grad_u(i, j) += dphi[M][QP](j) * Uib[i][M];

//  double epskk_ = 0;
//  for (uint k=0; k<3; k++ ) epskk_ += P->grad_u(k,k);

//  P->sigeff = 0;
//  for (uint i=0; i<3; i++) for (uint j=0; j<3; j++) 
//  for (uint k=0; k<3; k++) for (uint l=0; l<3; l++)
//    P->sigeff(i,j) += C_ijkl(i,j,k,l) * P->grad_u(k,l);

//  // \sig_tot = sig_eff - \alpha_d * T * \delta_ij
//  P->sigtot = P->sigeff;
//  for (uint k=0; k<3; k++ ) 
//    P->sigtot(k,k) -= P->alpha_d * T->temperature;

//  // dev_ij = sig_ij - 1/3 \delta_ij \sigma_kk
//  P->deviatoric = P->sigtot;
//  for (uint i=0; i<3; i++ )
//  for (uint k=0; k<3; k++ ) 
//    P->deviatoric(i,i) -= (1./3.) * P->sigtot(k,k);

//  double J2=0;
//  for (uint i=0; i<3; i++ )
//  for (uint j=0; j<3; j++ ) 
//    J2 += (1./2.) * P->deviatoric(i,j) * P->deviatoric(i,j);

//  P->von_mises = sqrt( 3 * J2 );
}

/**
 *
 */
void ViscoPlasticMaterial::project_stress( Elem & elem_ )
{
  SCOPELOG(10);
  /// Update the stresses in the interface
  reinit( *(system.current_local_solution), elem_ );
  do { update_ifc_qp(); } while ( next_qp() );

  ///  Project into the stress system
  stress_postproc.reinit( *elem );
  vp_ifc.reinit( elem->id(), qrule.n_points() );

  ViscoplasticIFC::PropsTranspose Pt( vp_ifc.by_qp );
  stress_postproc.project_tensor( Pt.sigtot, "sigtot" );
  stress_postproc.project_tensor( Pt.sigeff, "sigeff" );
  stress_postproc.project( Pt.von_mises, "von_mises" );
}


/**
 *
 *
 */
ad::VectorXreal ViscoPlasticMaterial::residual_qp( const ad::VectorXreal & /* _ad_Uib */ )
{
  const vector<Real> & JxW = fe->get_JxW();
  const vector<vector<RealGradient>> & dphi = fe->get_dphi();
  uint n_dofsv = dof_indices_var[0].size();

  for (uint i=0; i<3; i++)
  for (uint j=0; j<3; j++) 
  for (uint M=0;  M<n_dofsv;  M++) 
    grad_u(i, j) += dphi[M][QP](j) * Uib(i,M); /** Jacobian **/ // ****

  // ( \phi_j , Cijkl U_k,l ) 
  for (uint B=0;  B<n_dofsv;  B++)
  for (uint i=0; i<3; i++) 
  for (uint j=0; j<3; j++) 
  for (uint k=0; k<3; k++) 
  for (uint l=0; l<3; l++) 
    Fib(i,B) += JxW[QP] *  dphi[B][QP](j) * C_ijkl(i,j,k,l) * grad_u(k,l) ;

  // Temperature
  //       alpha_d = beta_d * K_bulk
  //       - ( \phi_i,i , \alpha_d T ) ==> term in the RHS.
  for (uint B=0;  B<n_dofsv;  B++)
  for (uint i=0;  i<3;  i++)
    Fib(i,B) -= JxW[QP] *  dphi[B][QP](i) * P->alpha_d  * T->temperature;

//  // Subtract the plastic strain as a body force
//  //
//  // - ( \phi_j , Cijkl \varepsilon^p_kl ) 
//  for (uint B=0;  B<n_dofsv;  B++)
//  for (uint i=0; i<3; i++) 
//  for (uint j=0; j<3; j++) 
//  for (uint k=0; k<3; k++) 
//  for (uint l=0; l<3; l++) 
//    Fib(i,B) -= JxW[QP] *  dphi[B][QP](j) * C_ijkl(i,j,k,l) * P->plastic_strain(k,l) ;

  return _ad_Fib;
}

ad::VectorXreal foo( const ad::VectorXreal & x )
{
  ad::VectorXreal out(2); 
  out << x[0], 5; 
  return out;
}

/**
 *    Adds the contribution of th quadrature point _qp_ to the element matrices Ke_var and vector.
 *
 *    Note: the solution vector must be parsed during reinit.
 */
void ViscoPlasticMaterial::residual_and_jacobian_qp ()
{

//  auto f = [&](const ad::VectorXreal & x) -> ad::VectorXreal
//  { return this->residual_qp(x);  };
//  ad::VectorXreal x;
//  ad::VectorXreal F;
//  double y[20];
//  Eigen::MatrixXd JAC = ad::jacobian( foo, ad::wrt(x), ad::at(x), F );

}

/**
 *     Builds the RHS of the element and assembles in the global _residual_ and _jacobian_.
 */
void ViscoPlasticMaterial::residual_and_jacobian (Elem & elem_, const NumericVector<Number> & soln, SparseMatrix<Number> * jacobian , NumericVector<Number> * residual )
{
  SCOPELOG(5);
  
  // Reinit the object
  reinit( soln, elem_ );

  // Build the element jacobian and residual for each quadrature point _qp_.
  do { residual_and_jacobian_qp(); } while ( next_qp() );

  /*
  *    Build the Re vector and assemble in the global residual vector
  */
//  uint n_dofsv = dof_indices_var[0].size();
//  for (uint i=0; i<3; i++) 
//  for (uint B=0;  B<n_dofsv;  B++)
//    Re( i*n_dofsv + B ) = Fib(i,B);

  const DofMap & dof_map = system.get_dof_map();

  if ( residual ) 
  {
//    dlog(1) << "Adding vector to residual ..." << Re;
    dof_map.constrain_element_vector ( Re, dof_indices );
    residual->add_vector ( Re, dof_indices );
  }

  // Add the the global matrix
  if ( jacobian ) 
  {
//    dlog(1) << "Adding matrix to jacobian ..." << Ke;
    dof_map.constrain_element_matrix (Ke, dof_indices);
    jacobian->add_matrix (Ke, dof_indices);
  }
}

/**
 *
 *   BOUNDARY CONDITIONS
 *
 */

/**
 *
 */
void ViscoPlasticMaterialBC::residual_and_jacobian_qp ()
{
  SCOPELOG(1);
  dlog(1) << "QP:" << QP;
  const vector<Real> & JxW = fe->get_JxW();
  const vector<vector<Real>> & phi = fe->get_phi();
  const vector<Point> & normals = fe->get_normals(); 
  uint n_dofsv = dof_indices_var[0].size();

  if ( sigtot )
  for (uint B=0; B<n_dofsv; B++)
  for (uint i=0; i<3; i++) 
  for (uint j=0; j<3; j++) 
    Fib(i,B) -= JxW[QP] * (*sigtot)(i,j) * normals[QP](j) * phi[B][QP];
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
  uint n_dofsv = dof_indices_var[0].size();

//  // Build the Re vector
//  for (uint i=0; i<3; i++) 
//  for (uint B=0;  B<n_dofsv;  B++)
//    Re( i*n_dofsv + B ) = Fib[i][B];

  if ( residual )
  {
    const DofMap & dof_map = system.get_dof_map();
    dof_map.constrain_element_vector (Re, dof_indices);
    residual->add_vector (Re, dof_indices);
  }
}

