
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
 *
 */
ViscoPlasticMaterial::ViscoPlasticMaterial( suint sid_,
                                            const MaterialConfig & config, 
                                            TransientNonlinearImplicitSystem & sys_,
                                            ViscoplasticSolver & vpsolver_ ) :
  Material( sid_, config ),
  dof_indices_var(3),
  P( 0 ), T(0), bc_material(0),
  vpsolver(vpsolver_), system( sys_ ), 
  stress_system( vpsolver.stress_system ),
  stress_postproc( *this, stress_system )
{
  setup_variables();
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
  SCOPELOG(1);
  // Calc xyz (fe must have already been initialized)
  const std::vector<Point> & xyz = fe->get_xyz();

  // We have it in cache? Return;
  if ( vp_ifc.valid )  {
    dlog(1) << "Skipping...";
    return; 
  }

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
  if ( ! bc_material ) 
  {
    bc_material = new ViscoPlasticMaterialBC( sid, config, system, vpsolver );
    bc_material->init_fem();
  }

  return bc_material;
}

/**
 *
 */
ViscoPlasticMaterialBC::ViscoPlasticMaterialBC( suint sid_, const MaterialConfig & config_,
                                                TransientNonlinearImplicitSystem & sys_,
                                                ViscoplasticSolver & vpsolver_ ) :
  ViscoPlasticMaterial( sid_, config_, sys_, vpsolver_ )
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
  SCOPELOG(1);

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
void ViscoPlasticMaterial::reinit( const NumericVector<Number> & soln, const Elem & elem, uint side )
{
  SCOPELOG(1);

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
}

/**
 *     Updates the interface based on the current Uib solution
 */
void ViscoPlasticMaterial::update_ifc_qp()
{
  uint n_dofsv = dof_indices_var[0].size();
  const vector<vector<RealGradient>> & dphi = fe->get_dphi();

  // Computes grad_u  ==> move to reinit, and to the interface (?)
  P->grad_u = 0;
  for (uint i=0; i<3; i++)
  for (uint j=0; j<3; j++)
  for (uint M=0;  M<n_dofsv;  M++)
    P->grad_u(i, j) += dphi[M][QP](j) * Uib[i][M];

  double epskk_ = 0;
  for (uint k=0; k<3; k++ ) epskk_ += P->grad_u(k,k);

  P->sigeff = 0;
  for (uint i=0; i<3; i++) for (uint j=0; j<3; j++) 
  for (uint k=0; k<3; k++) for (uint l=0; l<3; l++)
    P->sigeff(i,j) += C_ijkl(i,j,k,l) * P->grad_u(k,l);

  // \sig_tot = sig_eff - \alpha_d * T * \delta_ij
  P->sigtot = P->sigeff;
  for (uint k=0; k<3; k++ ) 
    P->sigtot(k,k) -= P->alpha_d * T->temperature;

  // dev_ij = sig_ij - 1/3 \delta_ij \sigma_kk
  P->deviatoric = P->sigtot;
  for (uint i=0; i<3; i++ )
  for (uint k=0; k<3; k++ ) 
    P->deviatoric(i,i) -= (1./3.) * P->sigtot(k,k);

  double J2=0;
  for (uint i=0; i<3; i++ )
  for (uint j=0; j<3; j++ ) 
    J2 += (1./2.) * P->deviatoric(i,j) * P->deviatoric(i,j);

  P->von_mises = sqrt( 3 * J2 );
}

/**
 *
 */
void ViscoPlasticMaterial::project_stress()
{
  stress_postproc.reinit( *elem );
  vp_ifc.reinit( elem->id(), qrule.n_points() );

  ViscoPlasticIFC::PropsTranspose Pt( vp_ifc.by_qp );
  stress_postproc.project_tensor( Pt.sigtot, "sigtot" );
  stress_postproc.project_tensor( Pt.sigeff, "sigeff" );
  stress_postproc.project( Pt.von_mises, "von_mises" );
}

///**
// *     Calculates the output information at a single point in this material.
// *     Pushes the results into the trg_coupler
// *     _elem_ is the element in this mesh where the point lies
// */
//void ViscoPlasticMaterial::feed_coupler( ElemCoupler & trg_ec, const Point & trg_pt,
//                                         const Elem * elem, const NumericVector<Number> & soln )
//{
//  SCOPELOG(5);

//  // Interpolators are the ones from this material.
//  // However, they must be initialized with the xyz of the gauss points of the target
//  const vector<vector<Real>> & phi = fe->get_phi();
//  const vector<vector<RealGradient>> & dphi = fe->get_dphi();
//  const DofMap & dof_map = system.get_dof_map();
//  dof_map.dof_indices (elem, dof_indices);
//  dof_map.dof_indices (elem, dof_indices_var[0], 0);
//  uint n_dofsv = dof_indices_var[0].size();

//  std::vector<Point> pts_from = { trg_pt };
//  std::vector<Point> pts_to;
//  FEInterface::inverse_map (3, fe->get_fe_type(), elem, pts_from, pts_to, 1E-10);
//  fe->reinit( elem, & pts_to );

//  // Prepare the Uib vector for the automatic differentiation
//  Uib.clear();
//  for ( uint i=0; i<3; i++ )
//  {
//    vector<Number> row;
//    for ( uint B=0; B<n_dofsv; B++ )
//      row.push_back( soln( dof_indices[i*n_dofsv + B] ) );
//    Uib.push_back(row);
//  }

//  RealVectorValue U = 0;
//  for ( uint B=0; B<n_dofsv; B++ )
//  for ( uint i=0; i<3; i++ )
//    U(i) += phi[B][0] * Uib[i][B];

//  RealTensor GRAD;
//  for ( uint B=0; B<n_dofsv; B++ )
//  for ( uint i=0; i<3; i++ )
//  for ( uint j=0; j<3; j++ )
//    GRAD(i,j) += dphi[B][0](j) * Uib[i][B];

//  // Feed the coupler in this point
//  vector<RealVectorValue> & trg_Uqi  = trg_ec.vector_params["U"];
//  trg_Uqi.push_back( U );

//  vector<RealTensor> & trg_GRAD_Uqij = trg_ec.tensor_params["GRAD_U"];
//  trg_GRAD_Uqij.push_back( GRAD );

//}

/**
 *
 */
//void ViscoPlasticMaterial::update_plastic_strain()
//{
//  SCOPELOG(5) ;

////  dlog(1) << "Update Plastic Strain data:";
////  dlog(1) << "    Von Mises:                 " << von_mises;
////  dlog(1) << "    Deviatoric:                " << deviatoric;
////  dlog(1) << "    Plastic Strain (t):        " << plastic_strain;
////  dlog(1) << "    Delta t:                   " << ts.dt;

////  dlog(1) << "    Creep (Carter):" << ts.dt;
////  dlog(1) << "          A:                   " << creep_carter_a;
////  dlog(1) << "          Q:                   " << creep_carter_q;
////  dlog(1) << "          N:                   " << creep_carter_n;

//  // Update the values in the coupler (not in this object!)
//  auto & ec_plastic_strain = elem_coupler->tensor_params["plastic_strain"];
//  auto & ec_plastic_strain_rate = elem_coupler->tensor_params["plastic_strain_rate"];


//  // Some validation
//  uint nqp = qrule.n_points();
//  if ( deviatoric.size() != nqp )      flog << "Vector sizes do not match [deviatoric.size(" << deviatoric.size() << ") != nqp(" << nqp << ")]";
//  if ( ec_plastic_strain.size() != nqp )  flog << "Vector sizes do not match [plastic_strain.size(" << ec_plastic_strain.size() << ") != nqp(" << nqp << ")]";
//  if ( creep_carter_a.size() != nqp )  flog << "Vector sizes do not match [creep_carter_a.size(" << creep_carter_a.size() << ") != nqp(" << nqp << ")]";
//  if ( creep_carter_q.size() != nqp )  flog << "Vector sizes do not match [creep_carter_q.size(" << creep_carter_q.size() << ") != nqp(" << nqp << ")]";
//  if ( creep_carter_n.size() != nqp )  flog << "Vector sizes do not match [creep_carter_n.size(" << creep_carter_n.size() << ") != nqp(" << nqp << ")]";
//  if ( von_mises.size() != nqp )       flog << "Vector sizes do not match [von_mises.size(" << von_mises.size() << ") != nqp(" << nqp << ")]";
//  if ( temperature.size() != nqp )     flog << "Vector sizes do not match [temperature.size(" << temperature.size() << ") != nqp(" << nqp << ")]";

//  double R_ = 8.3144;   // Universal gas constant [ J/mol/K ]
//  ec_plastic_strain_rate.clear();
//  for ( uint qp=0; qp<qrule.n_points() ; qp++ )
//  {
//    double A_ = P.creep_carter_a, Q_=P.creep_carter_q, N_=P.creep_carter_n;

//    RealTensor psr = 3./2. * A_ * exp( - Q_ / R_ / TProp.temperature ) * 
//                     pow( P.von_mises/1e6 , N_-1 ) *
//                     P.deviatoric/1e6;

////    dlog(1) << "COEF1:" << 3./2. * A_;
////    dlog(1) << "COEF2:" << Q_ / R_ / temperature[qp];
////    dlog(1) << "COEF3:" << exp( - Q_ / R_ / temperature[qp] );
////    dlog(1) << "TEMP:" <<   temperature[qp];
////    dlog(1) << "POW VM:" << pow( von_mises[qp]/1e6 , N_-1 );
////    dlog(1) << "DEVIATORIC (MPa):" << deviatoric[qp]/1e6;
////    dlog(1) << "PSR: " << psr;
////    ec_plastic_strain_rate.push_back( psr );
////    ec_plastic_strain[qp] += psr * ts.dt;
//  }
////  dlog(1) << "    Plastic Strain rate:        " << plastic_strain_rate;
////  dlog(1) << "    Plastic Strain:             " << ec_plastic_strain;
////  dlog(1) << "    dt :                        " << ts.dt;
//}


/**
 *    Adds the contribution of th quadrature point _qp_ to the element matrices Ke_var and vector.
 *
 *    Note: the solution vector must be parsed during reinit.
 */
void ViscoPlasticMaterial::residual_and_jacobian_qp ()
{
  dlog(1) << "QP=" << QP;
  const vector<Real> & JxW = fe->get_JxW();
  const vector<vector<RealGradient>> & dphi = fe->get_dphi();
  uint n_dofsv = dof_indices_var[0].size();

  // Computes grad_u  ==> this must be a AD variable. Needs to have the right type.
  RealTensor grad_u;
  for (uint i=0; i<3; i++)
  for (uint j=0; j<3; j++)
  for (uint M=0;  M<n_dofsv;  M++)
    grad_u(i, j) += dphi[M][QP](j) * Uib[i][M];

  /** Jacobian **/
  // ****
  // Effective mechanics
  //       ( \phi_i,j , C_ijkl u_k,l )   ok
  for (uint B=0; B<n_dofsv;  B++)
  for (uint M=0; M<n_dofsv;  M++)
  for (uint i=0; i<3; i++) 
  for (uint j=0; j<3; j++)
  for (uint k=0; k<3; k++) 
  for (uint l=0; l<3; l++) 
  {
    Ke_var[i][k](B,M) += JxW[QP] * C_ijkl(i,j,k,l) * dphi[M][QP](l) * dphi[B][QP](j);
  }

  /** RESIDUAL **/

  // ( \phi_j , Cijkl U_k,l ) 
  for (uint B=0;  B<n_dofsv;  B++)
  for (uint i=0; i<3; i++) 
  for (uint j=0; j<3; j++) 
  for (uint k=0; k<3; k++) 
  for (uint l=0; l<3; l++) 
    Fib[i][B] += JxW[QP] *  dphi[B][QP](j) * C_ijkl(i,j,k,l) * grad_u(k,l) ;

  // Subtract the plastic strain as a body force
  //
  // - ( \phi_j , Cijkl \varepsilon^p_kl ) 
  for (uint B=0;  B<n_dofsv;  B++)
  for (uint i=0; i<3; i++) 
  for (uint j=0; j<3; j++) 
  for (uint k=0; k<3; k++) 
  for (uint l=0; l<3; l++) 
    Fib[i][B] -= JxW[QP] *  dphi[B][QP](j) * C_ijkl(i,j,k,l) * P->plastic_strain(k,l) ;

//  dlog(1) << "VPM::residual plastic_strain: " << plastic_strain;

  // Temperature
  //       alpha_d = beta_d * K_bulk
  //       - ( \phi_i,i , \alpha_d T ) ==> term in the RHS.
  for (uint B=0;  B<n_dofsv;  B++)
  for (uint i=0;  i<3;  i++)
    Fib[i][B] -= JxW[QP] *  dphi[B][QP](i) * P->alpha_d  * T->temperature;

}

/**
 *     Builds the RHS of the element and assembles in the global _residual_ and _jacobian_.
 */
void ViscoPlasticMaterial::residual_and_jacobian (Elem & elem, const NumericVector<Number> & soln, SparseMatrix<Number> * jacobian , NumericVector<Number> * residual )
{
  SCOPELOG(1);
  
  // Reinit the object
  reinit( soln, elem );

  // Build the element jacobian and residual for each quadrature point _qp_.
  do { residual_and_jacobian_qp(); } while ( next_qp() );

  /*
  *    Build the Re vector and assemble in the global residual vector
  */
  uint n_dofsv = dof_indices_var[0].size();
  for (uint i=0; i<3; i++) 
  for (uint B=0;  B<n_dofsv;  B++)
    Re( i*n_dofsv + B ) = Fib[i][B];

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
  SCOPELOG(5);
  const vector<Real> & JxW = fe->get_JxW();
  const vector<vector<Real>> & phi = fe->get_phi();
  const vector<Point> & normals = fe->get_normals(); 
  uint n_dofsv = dof_indices_var[0].size();

  if ( sigtot )
  for (uint B=0; B<n_dofsv; B++)
  for (uint i=0; i<3; i++) 
  for (uint j=0; j<3; j++) 
    Fib[i][B] -= JxW[QP] * (*sigtot)(i,j) * normals[QP](j) * phi[B][QP];
}

/**
 *     Builds the RHS of the element and assembles in the global _residual_.
 */
void ViscoPlasticMaterialBC::residual_and_jacobian ( Elem & elem, uint side, 
                                                     const NumericVector<Number> & soln, 
                                                     SparseMatrix<Number> * jacobian ,
                                                     NumericVector<Number> * residual )
{
  UNUSED(jacobian);
  // Reinit the object
  reinit( soln, elem, side );

  // Build the element jacobian and residual for each quadrature point _qp_.
  do { residual_and_jacobian_qp(); } while ( next_qp() );
  uint n_dofsv = dof_indices_var[0].size();

  // Build the Re vector
  for (uint i=0; i<3; i++) 
  for (uint B=0;  B<n_dofsv;  B++)
    Re( i*n_dofsv + B ) = Fib[i][B];

  if ( residual )
  {
    const DofMap & dof_map = system.get_dof_map();
    dof_map.constrain_element_vector (Re, dof_indices);
    residual->add_vector (Re, dof_indices);
  }
}


/**
 *
 */
ostream& operator<<(ostream& os, const ViscoPlasticIFC & m)
{
  os << "VISCOPLASTIC INTERFACE:" << endl;

  for ( auto & [ eid, pvec ] : m.by_elem )
  {
    os << "    EID:" << eid;
    os << "          [" << endl;
    for ( auto & p : pvec ) {
      os << "             lame_mu:" << setw(15) << p.lame_mu << endl;
      os << "         lame_lambda:" << setw(15) << p.lame_lambda << endl;
    }
    os << "          ]" << endl;
  }

  return os;
}
