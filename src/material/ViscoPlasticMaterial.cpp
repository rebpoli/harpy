
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
#include "util/Stopwatch.h"
#include "harpy/Timestep.h"


/**
 *    BC_OBJECT is true when this is being called from the children constructor
 */
ViscoPlasticMaterial::ViscoPlasticMaterial( suint sid_,
                                            const MaterialConfig * config, 
                                            TransientNonlinearImplicitSystem & sys_,
                                            ViscoplasticSolver & vpsolver_,
                                            bool called_from_bc_constructor ) :
  Material( sid_, config ),
  dof_indices_var(3),
  P( 0 ), bc_material(0),
  vpsolver(vpsolver_), system( sys_ ), 
  stress_system( vpsolver.stress_system ),
  stress_postproc( this, stress_system ),
  dfile( "run/csv/plasticity-"+fmt_i(RANK)+string("-sid_") + fmt_i(sid) + string("-") + config->name + "-" + string(called_from_bc_constructor?"BC":"Body") + string(".csv") ),
  res_jac_k(0)
{
  // Setup cariables only if on the parent
  if (! called_from_bc_constructor ) 
    bc_material = new ViscoPlasticMaterialBC( sid, config, system, vpsolver );

  dfile << "sid" << "timestep" << "res_jac_k" << "dt" << "newton_k" << "elem" << "qp" << "von_mises";
  dfile << "eps_kk" << "gradu" << "sigeff" << "sigtot" << "deviatoric" << "sigeff0" << "plast0";
  dfile << "F" << "PlStrain_n" << "PlStrain" << "PlStrRate" << "eps_tr" << "eps_tr*" << "zeta" << "U";
  dfile << endrow;
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
  // Calc xyz (fe must have already been initialized)
  const std::vector<Point> & xyz = fe->get_xyz();

  // We have it in cache? Return;
  if ( vp_ifc.valid ) return; 

  uint qp = 0;
  for ( auto & pt : xyz )
  {
    auto & prop = vp_ifc.get( qp++ );
    prop.init_from_config( *config, pt );
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
ViscoPlasticMaterialBC::ViscoPlasticMaterialBC( suint sid_, const MaterialConfig * config_,
                                                TransientNonlinearImplicitSystem & sys_,
                                                ViscoplasticSolver & vpsolver_ ) :
  ViscoPlasticMaterial( sid_, config_, sys_, vpsolver_, true )
{
}


/**
 *  This can only be done after EquationSystems init.
 *  This is only done once for the material.
 */
void ViscoPlasticMaterial::init_fem()
{
  SCOPELOG(1);
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

  // Initialize the dependent materials/posprocs

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

  n_dofs = dof_indices.size();
  n_dofsv = dof_indices_var[0].size();

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
  // Initialize the flattened containers
  _init_autodiff();

  // Feed Uib with current solution
  for ( uint i=0; i<3; i++ )
  for ( uint B=0; B<n_dofsv; B++ )
    Uib(i,B) = soln( dof_indices[i*n_dofsv + B] );

}

/**
 *
 *  Rewind the plastic strain due to a timestep cut
 */
void ViscoPlasticMaterial::rewind( Elem & elem_ )
{
  reinit(elem_);
  do { 
    P->plastic_strain = P->plastic_strain_n; 
    P->plastic_strain_k = P->plastic_strain_n;
    P->creep_md.etr = P->creep_md.etr_n;
  } while ( next_qp() );

  // Update the probes accordingly   ////PROBE
  // This can go into the IFC interface...
  for ( auto & [ _, m1 ] : vp_ifc.probes_by_pname_by_elem )
  for ( auto & probe_ifc : m1[elem->id()] ) 
  {
    probe_ifc->props.plastic_strain = probe_ifc->props.plastic_strain_n;
    probe_ifc->props.creep_md.etr = probe_ifc->props.creep_md.etr_n;
  }
}

/**
 *
 */
void ViscoPlasticMaterial::project_stress( Elem & elem_ )
{
  SCOPELOG(10);

  // Update plastic strain history -- we could make it blindly in the interface,
  // without element awareness
  reinit(elem_);
  do { 
    P->plastic_strain_n = P->plastic_strain; 
    P->creep_md.etr_n = P->creep_md.etr;
  } while ( next_qp() );

  // Update the probes accordingly   ////PROBE
  // This can go into the IFC interface...
  for ( auto & [ _, m1 ] : vp_ifc.probes_by_pname_by_elem )
  for ( auto & probe_ifc : m1[elem->id()] ) 
  {
    probe_ifc->props.plastic_strain_n = probe_ifc->props.plastic_strain;
    probe_ifc->props.creep_md.etr_n = probe_ifc->props.creep_md.etr;
  }

  ///  Project into the stress system
  stress_postproc.reinit( elem_ );
  vp_ifc.reinit( elem_.id(), qrule.n_points() );

  PropsTranspose Pt( vp_ifc.by_qp );
  stress_postproc.project_tensor( Pt.sigtot, "sigtot" );
  stress_postproc.project_tensor_invariants( Pt.sigtot , "S" );

  stress_postproc.project_tensor_invariants( Pt.sigeff_terz , "Sterz" );
  
  stress_postproc.project_tensor( Pt.sigeff, "sigeff" );
  stress_postproc.project_tensor( Pt.deviatoric, "deviatoric" );
  stress_postproc.project_tensor( Pt.plastic_strain, "plastic_strain" );
  stress_postproc.project_tensor( Pt.plastic_strain_rate, "plastic_strain_rate" );
  stress_postproc.project_tensor( Pt.initial_stress, "initial_stress" );
  stress_postproc.project_tensor( Pt.initial_strain, "initial_strain" );
  stress_postproc.project( Pt.von_mises, "von_mises" );
  stress_postproc.project( Pt.epskk, "epskk" );
  stress_postproc.project( Pt.F, "F" );
}

/**
 * 
 * This is the actual calculation of the residual using the AutoDiff types.
 *
 */
AD::Vec ViscoPlasticMaterial::residual_qp( const AD::Vec & /* ad_Uib */ )
{
  if ( ! P->lame_mu ) flog << "Null lame_mu? Elem:" << elem->id() << "  Qp:" << QP;

  const vector<Real> & JxW = fe->get_JxW();
  const vector<vector<RealGradient>> & dphi = fe->get_dphi();
  const vector<vector<Real>> & phi = fe->get_phi();


//  bool deb = 0;
//  if ( ( elem->id() == 36330 ) && ( ! QP ) ) deb = 1;

  ad_Fib.setZero(); 

  AD::Mat grad_u(3,3);
  for (uint i=0; i<3; i++)
  for (uint j=0; j<3; j++) 
  for (uint M=0;  M<n_dofsv;  M++) 
    grad_u(i, j) += dphi[M][QP](j) * Uib(i,M); /** Jacobian **/ // ****

//  dlog(1) << "grad_u? " << grad_u;
//  dlog(1) << "initial_strain? " << P->initial_strain;

  // ( \phi_j , Cijkl U_k,l - Cijkl U^0_k,l ) 
  for (uint B=0;  B<n_dofsv;  B++)
  for (uint i=0; i<3; i++) 
  for (uint j=0; j<3; j++) 
  for (uint k=0; k<3; k++) 
  for (uint l=0; l<3; l++) 
    Fib(i,B) += JxW[QP] *  dphi[B][QP](j) * P->C_ijkl(i,j,k,l) * grad_u(k,l) ;// - P->initial_strain(k,l) ) ;

  // Init sigeff with the initial stuff
  for (uint B=0;  B<n_dofsv;  B++)
  for (uint i=0; i<3; i++ )
  for (uint j=0; j<3; j++ ) 
    Fib(i,B) += JxW[QP] * dphi[B][QP](j) * P->initial_stress(i,j);

  // - ( \phi , \rho g )
  for (uint B=0;  B<n_dofsv;  B++)
    Fib(2,B) -= JxW[QP] *  phi[B][QP] * P->density * GRAVITY_FORCE;   // Assuming Z increases with depth (Z axis points downward)

  // Temperature
  //       alpha_d = beta_d * K_bulk
  //       - ( \phi_i,i , \alpha_d T ) ==> term in the RHS.
  for (uint B=0;  B<n_dofsv;  B++)
  for (uint i=0;  i<3;  i++)
    Fib(i,B) -= JxW[QP] *  dphi[B][QP](i) * P->alpha_d  * ( P->temperature - P->initial_temperature );

  // Pressure
  //       - ( \phi_i,i , biot P ) ==> term in the RHS.
  for (uint B=0;  B<n_dofsv;  B++)
  for (uint i=0;  i<3;  i++)
    Fib(i,B) -= JxW[QP] *  dphi[B][QP](i) * P->biot  * ( P->pressure - P->initial_pressure );

  // For visualization and debugging
  AD::real epskk = 0;
  for (uint k=0; k<3; k++ ) epskk += grad_u(k,k);

  // Update stresses and plastic strain
  AD::Mat sigeff(3,3);

  P->sigeff0 = AD::norm(sigeff);
  P->plast0 = P->plastic_strain_k.norm();
  P->plast0_t = P->plastic_strain_k;

  for (uint i=0; i<3; i++) 
  for (uint j=0; j<3; j++) 
    sigeff(i,j) += P->initial_stress(i,j);

//  if ( deb ) ilog << "plast_k:    " << P->plastic_strain_k.norm();

  for (uint i=0; i<3; i++) for (uint j=0; j<3; j++) 
  for (uint k=0; k<3; k++) for (uint l=0; l<3; l++)
    sigeff(i,j) += P->C_ijkl(i,j,k,l) * ( grad_u(k,l) - P->plastic_strain_k(k,l) );

  AD::Mat sigtot = sigeff;
  for (uint k=0; k<3; k++ ) 
    sigtot(k,k) -= P->alpha_d * ( P->temperature - P->initial_temperature );
  for (uint k=0; k<3; k++ ) 
    sigtot(k,k) -= P->biot * ( P->pressure - P->initial_pressure );

  AD::Mat deviatoric = sigtot;
  for (uint i=0; i<3; i++ )
  for (uint k=0; k<3; k++ ) 
    deviatoric(i,i) -= (1./3.) * sigtot(k,k);

  AD::real J2=0;
  for (uint i=0; i<3; i++ )
  for (uint j=0; j<3; j++ ) 
    J2 += (1./2.) * deviatoric(i,j) * deviatoric(i,j);

  AD::real von_mises = sqrt( 3 * J2 );

  double R_ = 8.3144;   // Universal gas constant [ J/mol/K ]

  // Steadystate crep
  AD::real creep_ss_rate = 0;
  for ( auto & ss : P->creep_md.ss )
    creep_ss_rate += exp( - pow(ss.q/R_/P->temperature, ss.stretch) ) *
                     pow( von_mises/ss.sig0 , ss.n );

  // Transient creep
  AD::real F = 1;
  for ( auto & tr : P->creep_md.tr )
  {
    AD::real etr_star = exp( tr.c * P->temperature ) * 
                        pow( ( von_mises/tr.sig0 ), tr.m );

    AD::real zeta = 0;
    if ( abs(val(etr_star)) > 1e-20 ) zeta = 1 - P->creep_md.etr / etr_star;

    double alpha = tr.alpha_w;
    if ( zeta <= 0 ) alpha = 0;
    F = F * exp( alpha * zeta * zeta );
  }

  AD::real etr_rate = ( F - 1 ) * creep_ss_rate;

  P->creep_md.etr = P->creep_md.etr_n + val(etr_rate) * vpsolver.ts.dt;

  /// NOTE: This must match the calculations in the update function (ViscoplasticIFC)
  // Compute plastic strain rate
  AD::Mat plastic_strain_rate = AD::Mat::Zero(3,3);

  // If von_mises=0, the strain rate is zero
  if ( von_mises )
  for ( auto & ss : P->creep_md.ss )
    plastic_strain_rate += 
              3./2. * F * 
              exp( - pow(ss.q/R_/P->temperature, ss.stretch) ) *
              pow( von_mises/ss.sig0 , ss.n-1 ) *
              deviatoric / ss.sig0;

  // Update the plasteic strain
  AD::Mat plastic_strain = vpsolver.ts.dt * plastic_strain_rate;
  for (uint i=0; i<3; i++ )
  for (uint j=0; j<3; j++ ) 
    plastic_strain(i,j) += P->plastic_strain_n(i,j);

  /// Set the results into the quadrature point
  AD::dump( sigtot, P->sigtot );
  AD::dump( sigeff, P->sigeff );
  AD::dump( deviatoric, P->deviatoric );
  AD::dump( plastic_strain_rate, P->plastic_strain_rate );
  AD::dump( plastic_strain, P->plastic_strain);
  AD::dump( grad_u, P->GRAD_U);

  // Debug
  const vector<Point> & xyz = fe->get_xyz();

  P->von_mises = val(von_mises);
  P->epskk = val(epskk);
  P->F = val(F);
  double g=0; 
  for (uint i=0; i<3; i++ ) for (uint j=0; j<3; j++ ) g += pow( val(grad_u(i,j)), 2 );
  P->grad_norm = g;

  // 
  // Subtract the plastic strain as a body force
  //
  // - ( \phi_j , Cijkl \varepsilon^p_kl ) 
  for (uint B=0;  B<n_dofsv;  B++)
  for (uint i=0; i<3; i++) 
  for (uint j=0; j<3; j++) 
  for (uint k=0; k<3; k++) 
  for (uint l=0; l<3; l++) 
    Fib(i,B) -= JxW[QP] * dphi[B][QP](j) * P->C_ijkl(i,j,k,l) * plastic_strain(k,l) ;

  //   Hunting NaN
  for (uint i=0; i<3; i++) 
  for (uint B=0;  B<n_dofsv;  B++)
    if ( P->plastic_strain.norm() > 1e30 ) 
    {
      ilog << "ELEM:" << elem->id();
      ilog << "QP:" << QP;
      ilog << "grad_u:" << Print(grad_u);
      ilog << "sigeff:" << Print(sigeff);
      ilog << "sigtot:" << Print(sigtot);
      ilog << "deviatoric:" << Print(deviatoric);
      ilog << "J2:" << J2;
      ilog << "K:" << P->creep_md;
      ilog << "T:" << P->temperature;
      ilog << "VM:" << von_mises;
      flog << "Found NaN in residual! " << *P;
    }
 
  return ad_Fib;
}

/**
 *    Adds the contribution of th quadrature point _qp_ to the element matrices Ke_var and vector.
 *
 *    Note: the solution vector must be parsed during reinit.
 */
void ViscoPlasticMaterial::residual_and_jacobian_qp ()
{
  // Lambda function to compatibilize stuff
  auto f = [this](const AD::Vec & x) { return this->residual_qp(x);  };

  AD::Vec F;
  ad_Jijbm = AD::jacobian( f, wrt(ad_Uib), at(ad_Uib), F );

  // Update plastic_strain_k : the plastic strain after this newton iteration
  P->plastic_strain_k = P->plastic_strain;

  dfile << sid;
  dfile << vpsolver.ts.t_step;
  dfile << res_jac_k;
  dfile << vpsolver.ts.dt;
  dfile << 0;
  dfile << elem->id();
  dfile << QP;
  dfile << P->von_mises;
  dfile << P->epskk;
  dfile << P->grad_norm;
  dfile << P->sigeff.norm();
  dfile << P->sigtot.norm();
  dfile << P->deviatoric.norm();
  dfile << P->sigeff0;
  dfile << P->plast0;
  dfile << P->F;
  dfile << P->plastic_strain_n.norm();
  dfile << P->plastic_strain.norm();
  dfile << P->plastic_strain_rate.norm();
  dfile << P->creep_md.etr;
  //Compute the sum of U
  double u = 0; 
  for ( uint i=0; i<3; i++ )
  for ( uint B=0; B<n_dofsv; B++ )
      u += pow( val(Uib(i,B)) , 2 ) ;
  dfile << u;
//  dfile << Print(P->plastic_strain);
//  dfile << Print(P->plast0_t);
  dfile << endrow;

  // Map from the AD variable to libmesh datastructures
  for (uint i=0; i<3; i++) for (uint B=0;  B<n_dofsv;  B++)
    if ( isnan( val(Fib(i,B) ) ) ) 
      flog << "Found NaN in residual! " << *P;

  // Search for nan
  // Map from the AD variable to libmesh datastructures
  for (uint B=0; B<n_dofsv;  B++) for (uint M=0; M<n_dofsv;  M++)
  for (uint i=0; i<3; i++) for (uint j=0; j<3; j++)
    if ( isnan( val(Jijbm(i,j,B,M) ) ) ) 
      flog << "Found NaN in jacobian!" << *P;

  // Map from the AD variable to libmesh datastructures
  for (uint B=0; B<n_dofsv;  B++)
  for (uint M=0; M<n_dofsv;  M++)
  for (uint i=0; i<3; i++) 
  for (uint j=0; j<3; j++)
    Ke_var[i][j](B,M) += val( Jijbm(i,j,B,M) );

  // Map from the AD variable to libmesh datastructures
  for (uint i=0; i<3; i++) 
  for (uint B=0;  B<n_dofsv;  B++)
    Re( i*n_dofsv + B ) += val( Fib(i,B) );
}




/**
 *     Builds the RHS of the element and assembles in the global _residual_ and _jacobian_.
 */
void ViscoPlasticMaterial::residual_and_jacobian (Elem & elem_, const NumericVector<Number> & soln, SparseMatrix<Number> * jacobian , NumericVector<Number> * residual )
{
  SCOPELOG(10);

  res_jac_k++; // Debugging stuff
  
  // Reinit the object
  reinit( soln, elem_ );

  // Build the element jacobian and residual for each quadrature point _qp_.
  do { residual_and_jacobian_qp(); } while ( next_qp() );

  // Add to the residual 
  if ( residual ) 
  {
    const DofMap & dof_map = system.get_dof_map();
    dof_map.constrain_element_vector ( Re, dof_indices );
    residual->add_vector ( Re, dof_indices );
  }

  // Add the the global matrix
  if ( jacobian ) 
  {
    const DofMap & dof_map = system.get_dof_map();
    dof_map.constrain_element_matrix (Ke, dof_indices);
    jacobian->add_matrix (Ke, dof_indices);
  }

  // Keep the probes in sync with the solution
  update_probes();  ///PROBE
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

  ad_Fib.setZero(); 

  // Note: as the sigtot is constant, we dont need the jacobian
  //     : Keeping the variable as a AUTODIFF for parallelism
  //     : with the parent object. 
  if ( sigtot )
  for (uint B=0; B<n_dofsv; B++)
  for (uint i=0; i<3; i++) 
  for (uint j=0; j<3; j++) 
    Fib(i,B) -= JxW[QP] * (*sigtot)(i,j) * normals[QP](j) * phi[B][QP];

  // Map from the AD variable to libmesh datastructures
  for (uint i=0; i<3; i++) 
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

/**
 *
 */
void ViscoPlasticMaterial::apply_strain_initialization_method()
{
  using enum MatInitializeMethod;

  /** Hydrostatic ? **/
  if ( config->initialize.method == HYDROSTATIC ) 
  {
    ilog1 << "Forcing initial strain to hydrostatic in material '" << name << "'.";
    for ( auto & [ eid, pvec ] : vp_ifc.by_elem )
    for ( VPProps & p : pvec ) 
      p.force_initial_strain_to_hydrostatic();

    // Probes
    for ( auto & [ _, m1 ] : vp_ifc.probes_by_pname_by_elem )
    for ( auto & probe_ifc : m1[elem->id()] ) 
      probe_ifc->props.force_initial_strain_to_hydrostatic();
  }

}

/**
 *  PUT THE INITIAL STRAIN IN PLACE 
 */
void ViscoPlasticMaterial::update_initial_strain()
{
  // Gauss points of the elements
  for ( auto & [ eid, pvec ] : vp_ifc.by_elem )
  for ( VPProps & p : pvec ) 
    p.update_initial_strain();

  // Probes
  for ( auto & [ _, m1 ] : vp_ifc.probes_by_pname_by_elem )
  for ( auto & probe_ifc : m1[elem->id()] ) 
    probe_ifc->props.update_initial_strain();
  
  /** Do the necessary adjustments **/
  apply_strain_initialization_method();
}

/**
 *  This should be called before the probes export data.
 *
 */
void ViscoPlasticMaterial::update_probes()
{
  for ( auto & [ _, m1 ] : vp_ifc.probes_by_pname_by_elem )
  for ( auto & probe_ifc : m1[elem->id()] ) 
  {
    const Point & pt = probe_ifc->pt;
    auto & props = probe_ifc->props;
    props_at( props, pt, elem );
  }
}

/**
 *     Calculates the output information at a single point in this material.
 *     Pushes the results into the trg_coupler
 *     _elem_ is the element in this mesh where the point lies
 */
void ViscoPlasticMaterial::props_at( VPProps & p, 
                                     const Point & pt, const Elem * elem )
{
  SCOPELOG(5);
  const DofMap & dof_map = system.get_dof_map();
  const vector<vector<Real>> & phi = fe->get_phi();
  const vector<vector<RealGradient>> & dphi = fe->get_dphi();

  // Calculate Uib (avoiding reusage to reduce coupling)
  dof_map.dof_indices (elem, dof_indices_var[0], 0);
  dof_map.dof_indices (elem, dof_indices_var[1], 1);
  dof_map.dof_indices (elem, dof_indices_var[2], 2);
  n_dofsv = dof_indices_var[0].size();
  vector< vector<double> > _Uib(3, vector<double>(n_dofsv,0));

  // Reinit FE in the point we want
  std::vector<Point> pts_from = { pt };
  std::vector<Point> pts_to;
  FEInterface::inverse_map (3, fe->get_fe_type(), elem, pts_from, pts_to, 1E-10);
  fe->reinit( elem, & pts_to );

  // Compute U and GRAD_U at the point
  RealVectorValue U = 0;
  for ( uint B=0; B<n_dofsv; B++ )
  for ( uint i=0; i<3; i++ )
    U(i) += phi[B][0] * val( Uib(i,B) );

  RealTensor GRAD_U;
  for ( uint B=0; B<n_dofsv; B++ )
  for ( uint i=0; i<3; i++ )
  for ( uint j=0; j<3; j++ )
    GRAD_U(i,j) += dphi[B][0](j) * val( Uib(i,B) );

  // Update the VPProps stresses, plasticity etc
  p.update( U, GRAD_U, vpsolver.ts.dt );
}



