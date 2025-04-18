
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
#include "harpy/Timestep.h"


/**
 *
 */
ViscoPlasticMaterial::ViscoPlasticMaterial( suint sid_,
                                            const MaterialConfig & config, 
                                            TransientNonlinearImplicitSystem & sys_,
                                            const Timestep & ts_) :
  Material( sid_, config ), dof_indices_var(3),
  system( sys_ ), implicit(1), bc_material(0), ts(ts_)
{
  // Lists the necessary properties to fetch from the config during init_coupler
  required_material_properties.assign({
      "porothermoelastic.alpha_d",
      "porothermoelastic.beta_e",
      "porothermoelastic.young",
      "porothermoelastic.poisson",
      "porothermoelastic.lame_mu",
      "porothermoelastic.lame_lambda",
      "creep_carter.a",
      "creep_carter.q",
      "creep_carter.n"
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
    bc_material = new ViscoPlasticMaterialBC( sid, config, system, ts );
    bc_material->init_fem();
  }

  return bc_material;
}

/**
 *
 */
ViscoPlasticMaterialBC::ViscoPlasticMaterialBC( suint sid_, const MaterialConfig & config_,
                                                TransientNonlinearImplicitSystem & sys_,
                                                const Timestep & ts_) :
  ViscoPlasticMaterial( sid_, config_, sys_, ts_ )
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
}

/**
 *  Init the DoF map and the element matrices.
 *
 */
void ViscoPlasticMaterial::reinit( Coupler & coupler, const Elem & elem, uint side )
{
  SCOPELOG(3);

  /*
   * Initializes the FEM structures, not depending on the solution
   */
  reinit( elem, side );

  elem_coupler = & ( coupler.elem_coupler(elem.id()) );

  fetch_from_coupler();
}

/**
 *
 */
void ViscoPlasticMaterial::fetch_from_coupler()
{
  SCOPELOG(5);

  // Fetch the needed parameters from the coupler
  // From the configuration file
  get_from_element_coupler(  "porothermoelastic.alpha_d",       alpha_d        ); 
  get_from_element_coupler(  "porothermoelastic.lame_mu",       lame_mu        ); 
  get_from_element_coupler(  "porothermoelastic.lame_lambda",   lame_lambda    ); 

  // Creep: Carter model
  get_from_element_coupler(  "creep_carter.a",   creep_carter_a    ); 
  get_from_element_coupler(  "creep_carter.q",   creep_carter_q    ); 
  get_from_element_coupler(  "creep_carter.n",   creep_carter_n    ); 

  // From the thermal solver
  get_from_element_coupler(  "T",             temperature    ); 

  // From the stress solver
  get_from_element_coupler(  "von_mises",     von_mises      ); 
  get_from_element_coupler(  "deviatoric",    deviatoric     ); 

  // Direct access to makes sure our plastic_strain is initialized in the coupler 
  // Makes sure it has zeros for each qp if not initialized
  auto & ec_plastic_strain = elem_coupler->tensor_params["plastic_strain"]; 
  if  ( ! ec_plastic_strain.size() ) 
  for ( uint qp=0; qp<qrule.n_points() ; qp++ )
    ec_plastic_strain.push_back(RealTensor());

  // Direct access to makes sure our plastic_strain is initialized in the coupler 
  // This is only an intermediate calculation, important for debugging and visualization only
  plastic_strain = elem_coupler->tensor_params["plastic_strain"]; 
  plastic_strain_rate = elem_coupler->tensor_params["plastic_strain_rate"]; 
}

/**
 *  Init the DoF map and the element matrices.
 *
 *  Then, initializes the solution depending on soln.
 */
void ViscoPlasticMaterial::reinit( const NumericVector<Number> & soln, Coupler & coupler, const Elem & elem, uint side )
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

  fetch_from_coupler();
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
 *
 */
void ViscoPlasticMaterial::update_plastic_strain()
{
  SCOPELOG(5) ;

//  dlog(1) << "Update Plastic Strain data:";
//  dlog(1) << "    Von Mises:                 " << von_mises;
//  dlog(1) << "    Deviatoric:                " << deviatoric;
//  dlog(1) << "    Plastic Strain (t):        " << plastic_strain;
//  dlog(1) << "    Delta t:                   " << ts.dt;

//  dlog(1) << "    Creep (Carter):" << ts.dt;
//  dlog(1) << "          A:                   " << creep_carter_a;
//  dlog(1) << "          Q:                   " << creep_carter_q;
//  dlog(1) << "          N:                   " << creep_carter_n;

  // Update the values in the coupler (not in this object!)
  auto & ec_plastic_strain = elem_coupler->tensor_params["plastic_strain"];
  auto & ec_plastic_strain_rate = elem_coupler->tensor_params["plastic_strain_rate"];


  // Some validation
  uint nqp = qrule.n_points();
  if ( deviatoric.size() != nqp )      flog << "Vector sizes do not match [deviatoric.size(" << deviatoric.size() << ") != nqp(" << nqp << ")]";
  if ( ec_plastic_strain.size() != nqp )  flog << "Vector sizes do not match [plastic_strain.size(" << ec_plastic_strain.size() << ") != nqp(" << nqp << ")]";
  if ( creep_carter_a.size() != nqp )  flog << "Vector sizes do not match [creep_carter_a.size(" << creep_carter_a.size() << ") != nqp(" << nqp << ")]";
  if ( creep_carter_q.size() != nqp )  flog << "Vector sizes do not match [creep_carter_q.size(" << creep_carter_q.size() << ") != nqp(" << nqp << ")]";
  if ( creep_carter_n.size() != nqp )  flog << "Vector sizes do not match [creep_carter_n.size(" << creep_carter_n.size() << ") != nqp(" << nqp << ")]";
  if ( von_mises.size() != nqp )       flog << "Vector sizes do not match [von_mises.size(" << von_mises.size() << ") != nqp(" << nqp << ")]";
  if ( temperature.size() != nqp )     flog << "Vector sizes do not match [temperature.size(" << temperature.size() << ") != nqp(" << nqp << ")]";

  double R_ = 8.3144;   // Universal gas constant [ J/mol/K ]
  ec_plastic_strain_rate.clear();
  for ( uint qp=0; qp<qrule.n_points() ; qp++ )
  {
    double A_ = creep_carter_a[qp], Q_=creep_carter_q[qp], N_=creep_carter_n[qp];

    RealTensor psr = 3./2. * A_ * exp( - Q_ / R_ / temperature[qp] ) * 
                     pow( von_mises[qp]/1e6 , N_-1 ) *
                     deviatoric[qp]/1e6;

//    dlog(1) << "COEF1:" << 3./2. * A_;
//    dlog(1) << "COEF2:" << Q_ / R_ / temperature[qp];
//    dlog(1) << "COEF3:" << exp( - Q_ / R_ / temperature[qp] );
//    dlog(1) << "TEMP:" <<   temperature[qp];
//    dlog(1) << "POW VM:" << pow( von_mises[qp]/1e6 , N_-1 );
//    dlog(1) << "DEVIATORIC (MPa):" << deviatoric[qp]/1e6;
//    dlog(1) << "PSR: " << psr;
    ec_plastic_strain_rate.push_back( psr );
    ec_plastic_strain[qp] += psr * ts.dt;
  }
//  dlog(1) << "    Plastic Strain rate:        " << plastic_strain_rate;
//  dlog(1) << "    Plastic Strain:             " << ec_plastic_strain;
//  dlog(1) << "    dt :                        " << ts.dt;
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

  // Subtract the plastic strain as a body force
  //
  // - ( \phi_j , Cijkl \varepsilon^p_kl ) 
  for (uint qp=0; qp<qrule.n_points(); qp++   )
  for (uint B=0;  B<n_dofsv;  B++)
  for (uint i=0; i<3; i++) 
  for (uint j=0; j<3; j++) 
  for (uint k=0; k<3; k++) 
  for (uint l=0; l<3; l++) 
    Fib[i][B] -= JxW[qp] *  dphi[B][qp](j) * C_ijkl(qp, i,j,k,l) * plastic_strain[qp](k,l) ;

//  dlog(1) << "VPM::residual plastic_strain: " << plastic_strain;

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


