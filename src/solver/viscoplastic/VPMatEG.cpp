#include "solver/viscoplastic/VPMatEG.h"

#include "config/ModelConfig.h"
#include "config/MaterialConfig.h"
#include "solver/viscoplastic/VPSolver.h"
#include "util/OutputOperators.h"

#include "libmesh/mesh.h"
#include "libmesh/dof_map.h"
#include "libmesh/fe_map.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/dense_submatrix.h"
#include "libmesh/dense_subvector.h"

namespace solver {
namespace viscoplastic {

VPMatEG::VPMatEG( TransientNonlinearImplicitSystem & sys_,
                  ViscoplasticSolver & vpsolver_ ) :
system(sys_),
vpsolver(vpsolver_),
fem_p(system),
fem_n(system),
qrule(2),
Pp(0), Pn(0),
n_uvars( vpsolver.is_eg() ? 6 : 3 ),
n_dofs_cg(0), n_dofs_eg(0),
QP(0),
ad(n_uvars)
{ 
  ASSERT( vpsolver.is_eg() , "Must be an EG solver to be here." );
  ASSERT( system.has_variable("UegX"), "Must be an EG solver to be here." );

  // Init Qrule
  uint vid = system.variable_number( "UegX" );
  DofMap & dof_map = system.get_dof_map();
  FEType fe_type = dof_map.variable_type(vid);
  qrule = QGauss( 2, fe_type.default_quadrature_order() );

  // Do not attach quadrature in fem_n -- it should be mapped by 
  // points during reinit (slave mapping)
  fem_p.attach_qrule( & qrule );
}

/**
 *
 */
void VPMatEG::init_properties()
{
  SCOPELOG(1);
  for ( EGFacePair & fp : gamma_I )
  {
    const MaterialConfig & cfg_p = vpsolver.get_material_config( fp.eid_p );
    const MaterialConfig & cfg_n = vpsolver.get_material_config( fp.eid_n );

    reinit( fp );

    const std::vector<Point> & xyz_p = fem_p.fe->get_xyz();
    const std::vector<Point> & xyz_n = fem_n.fe->get_xyz();
    for ( uint qp=0 ; qp<xyz_n.size() ; qp++ )
    {
      ASSERT( (xyz_p[qp]-xyz_n[qp]).norm() < 1e-8 , "These points should coincide.");
      const Point & pt = xyz_p[qp];
      fp.get_Pp(qp)->init_from_config( cfg_p, pt );
      fp.get_Pn(qp)->init_from_config( cfg_n, pt );
    }

    dlog(1) << fp;
  }
}

/**
 *
 */
void VPMatEG::init()
{
  SCOPELOG(1);
  MeshBase & mesh = system.get_mesh();

  for ( const auto & elem_p : mesh.active_local_element_ptr_range() )
  for ( auto side_p : elem_p->side_index_range() ) 
  {
    if ( elem_p->neighbor_ptr(side_p) != nullptr )
    {   // Register internal skeleton
      const Elem * elem_n = elem_p->neighbor_ptr(side_p);

      // Condition to add the penalty (avoid duplicates)
      //      Higher refinement level: prefer to integrate the lower levels
      //      Same refinement level: get the smaller id
      if ( ! elem_n->active() ) continue;                        // Inactive neighbor, continue.
      else if ( elem_p->level() > elem_n->level() ) continue;    // Higher level P. Continue
        else if (   ( elem_n->level() == elem_p->level() )       // Same level ?
                 && ( elem_p->id() > elem_n->id() ) ) continue;  // Lower ID wins.
      
      EGFacePair egfp( elem_p->id(), side_p, elem_n->id() );

      fem_p.fe->reinit( elem_p, side_p );
      uint nqp = qrule.n_points();

      egfp.Pq_p.resize(nqp);
      egfp.Pq_n.resize(nqp);

      gamma_I.emplace_back( egfp );
    }
    else
    {   // Register outer boundary skeleton
      EGFace egf( elem_p->id(), side_p );

      uint nqp = qrule.n_points();
      egf.Pq.resize(nqp);

      gamma_H.emplace_back( egf );
    }
  }

  init_properties();
}

/* Reinit structures for the face pair fp */
void VPMatEG::reinit( EGFacePair & fp )
{
  MeshBase & mesh = system.get_mesh();
  Elem * elem_p = mesh.elem_ptr(fp.eid_p);
  Elem * elem_n = mesh.elem_ptr(fp.eid_n);

  /* 1. Init FEM structs */
  const std::vector<Point> & xyz_p = fem_p.fe->get_xyz();
  fem_p.fe->reinit( elem_p , fp.side_p );

  vector<Point> pts;
  FEMap::inverse_map (3, elem_n, xyz_p, pts, 1E-10);
  fem_n.fe->reinit(elem_n, & pts);
}

/* Reinit structures for the face pair fp */
void VPMatEG::reinit( const NumericVector<Number> & soln , EGFacePair & fp )
{
  /* 1. Basic reinit */
  reinit(fp);

  /* 2. Update dofmaps and dof counters. */
  setup_dofs( fp );

  uint nd = dof_indices_eg.size();
  Ke.resize (nd, nd);
  Re.resize (nd);

  /* 3. Init QP counter and resolve autodiff solution vectors */
  QP = 0;

  // TODO: at some point the elements might have different number of dofs 
  //     : Only EG dofs
  ad.init( 0, n_dofs_eg, 2 ); 

  /* 4. Feed Uib with current solution - just EG part */
  for ( uint e=0; e<2; e++ )
  for ( uint i=3; i<6; i++ )
  for ( uint B=0; B<n_dofs_eg; B++ )
    ad.Ueib(e,i,B) = soln( dof_indices_eg[ad.idx(e,i,B)] );

  /* 5. Feed Uib with current solution - CG part */
  Ucg_eib.reserve( idx_cg(2,3,n_dofs_cg)+1 );
  Ucg_eib.clear();
  for ( uint e=0; e<2; e++ )
  for ( uint i=0; i<3; i++ )
  for ( uint B=0; B<n_dofs_cg; B++ )
  {
    uint ii = idx_cg(e,i,B);
    uint dof = dof_indices_cg[ii];
    Ucg_eib[ii] = soln(dof);
  }
}

/**
 *
 */
void VPMatEG::setup_dofs( EGFacePair & fp )
{
  MeshBase & mesh = system.get_mesh();
  Elem * elem_p = mesh.elem_ptr(fp.eid_p);
  Elem * elem_n = mesh.elem_ptr(fp.eid_n);

  /* Feed the dof indices in our datastructures */
  fem_p.set_dofs( system, elem_p );
  fem_n.set_dofs( system, elem_n );

  // Flat structure with EG dofs - for the dof assemblage
  dof_indices_eg.clear();
  dof_indices_eg.reserve( fem_p.dofi_eg.size() + fem_n.dofi_eg.size() );
  dof_indices_eg.insert( dof_indices_eg.end(), fem_p.dofi_eg.begin() , fem_p.dofi_eg.end() );
  dof_indices_eg.insert( dof_indices_eg.end(), fem_n.dofi_eg.begin() , fem_n.dofi_eg.end() );

  // Flat structure with CG dofs - for the dof assemblage
  dof_indices_cg.clear();
  dof_indices_cg.reserve( fem_p.dofi_cg.size() + fem_n.dofi_cg.size() );
  dof_indices_cg.insert( dof_indices_cg.end(), fem_p.dofi_cg.begin() , fem_p.dofi_cg.end() );
  dof_indices_cg.insert( dof_indices_cg.end(), fem_n.dofi_cg.begin() , fem_n.dofi_cg.end() );

  // Init dofs counters
  const DofMap & dof_map = system.get_dof_map();
  vector<dof_id_type> di;
  dof_map.dof_indices ( elem_p, di, 0 );
  n_dofs_cg = di.size();
  dof_map.dof_indices ( elem_p, di, 3 );
  n_dofs_eg = di.size();

  /* Useful validation */
  ASSERT( fp.Pq_p.size() == fp.Pq_n.size()    , "Size of the properties vector should be equal (" <<  fp.Pq_p.size() << " != " << fp.Pq_n.size() << ")" );
  ASSERT( fp.Pq_p.size() == qrule.n_points() , "Size of the properties vector should be equal to nqp (" << fp.Pq_p.size() << " != " << qrule.n_points() << ")" );
}

/**
 * 
 * This is the actual calculation of the residual using the AutoDiff types.
 *
 */
AD::Vec VPMatEG::residual_qp( const AD::Vec & /* ad_Uib */ )
{
  ASSERT( Pp->lame_mu || Pn->lame_mu, "Null lame_mu?" );
  return ad.ad_Fib;

//  const vector<Real> & JxW = fe->get_JxW();
//  const vector<vector<RealGradient>> & dphi = fe->get_dphi();
//  const vector<vector<Real>> & phi = fe->get_phi();


////  bool deb = 0;
////  if ( ( elem->id() == 36330 ) && ( ! QP ) ) deb = 1;

//  ad_Fib.setZero(); 

//  AD::Mat grad_u(3,3);
//  for (uint i=0; i<3; i++)
//  for (uint j=0; j<3; j++) 
//  for (uint M=0;  M<n_dofsv;  M++) 
//    grad_u(i, j) += dphi[M][QP](j) * Uib(i,M); /** Jacobian **/ // ****

//  // EG part - n_uvars > 3 only if EG is in use ; we could also use vpsolver.is_eg()
//  for (uint i=3; i<n_uvars; i++)
//  for (uint j=3; j<n_uvars; j++) 
//  for (uint M=0;  M<n_dofs_eg;  M++) 
//    grad_u(i, j) += dphi[M][QP](j) * Uib(i,M); /** Jacobian **/ // ****

//  // ( \phi_j , Cijkl U_k,l - Cijkl U^0_k,l ) 
//  for (uint B=0;  B<n_dofsv;  B++)
//  for (uint i=0; i<3; i++) 
//  for (uint j=0; j<3; j++) 
//  for (uint k=0; k<3; k++) 
//  for (uint l=0; l<3; l++) 
//    Fib(i,B) += JxW[QP] *  dphi[B][QP](j) * P->C_ijkl(i,j,k,l) * grad_u(k,l) ;

//  // Init sigeff with the initial stuff
//  for (uint B=0;  B<n_dofsv;  B++)
//  for (uint i=0; i<3; i++ )
//  for (uint j=0; j<3; j++ ) 
//    Fib(i,B) += JxW[QP] * dphi[B][QP](j) * P->initial_stress(i,j);

//  // - ( \phi , \rho g )
//  for (uint B=0;  B<n_dofsv;  B++)
//    Fib(2,B) -= JxW[QP] *  phi[B][QP] * P->density * GRAVITY_FORCE;   // Assuming Z increases with depth (Z axis points downward)

//  // Temperature
//  //       alpha_d = beta_d * K_bulk
//  //       - ( \phi_i,i , \alpha_d T ) ==> term in the RHS.
//  for (uint B=0;  B<n_dofsv;  B++)
//  for (uint i=0;  i<3;  i++)
//    Fib(i,B) -= JxW[QP] *  dphi[B][QP](i) * P->alpha_d  * ( P->temperature - P->initial_temperature );

//  // Pressure
//  //       - ( \phi_i,i , biot P ) ==> term in the RHS.
//  for (uint B=0;  B<n_dofsv;  B++)
//  for (uint i=0;  i<3;  i++)
//    Fib(i,B) -= JxW[QP] *  dphi[B][QP](i) * P->biot  * ( P->pressure - P->initial_pressure );

//  // For visualization and debugging
//  AD::real epskk = 0;
//  for (uint k=0; k<3; k++ ) epskk += grad_u(k,k);

//  // Update stresses and plastic strain
//  AD::Mat sigeff(3,3);

//  P->sigeff0 = AD::norm(sigeff);
//  P->plast0 = P->plastic_strain_k.norm();
//  P->plast0_t = P->plastic_strain_k;

//  for (uint i=0; i<3; i++) 
//  for (uint j=0; j<3; j++) 
//    sigeff(i,j) += P->initial_stress(i,j);

////  if ( deb ) ilog << "plast_k:    " << P->plastic_strain_k.norm();

//  for (uint i=0; i<3; i++) for (uint j=0; j<3; j++) 
//  for (uint k=0; k<3; k++) for (uint l=0; l<3; l++)
//    sigeff(i,j) += P->C_ijkl(i,j,k,l) * ( grad_u(k,l) - P->plastic_strain_k(k,l) );

//  AD::Mat sigtot = sigeff;
//  for (uint k=0; k<3; k++ ) 
//    sigtot(k,k) -= P->alpha_d * ( P->temperature - P->initial_temperature );
//  for (uint k=0; k<3; k++ ) 
//    sigtot(k,k) -= P->biot * ( P->pressure - P->initial_pressure );

//  AD::Mat deviatoric = sigtot;
//  for (uint i=0; i<3; i++ )
//  for (uint k=0; k<3; k++ ) 
//    deviatoric(i,i) -= (1./3.) * sigtot(k,k);

//  AD::real J2=0;
//  for (uint i=0; i<3; i++ )
//  for (uint j=0; j<3; j++ ) 
//    J2 += (1./2.) * deviatoric(i,j) * deviatoric(i,j);

//  AD::real von_mises = sqrt( 3 * J2 );

//  double R_ = 8.3144;   // Universal gas constant [ J/mol/K ]

//  // Steadystate crep
//  AD::real creep_ss_rate = 0;
//  for ( auto & ss : P->creep_md.ss )
//    creep_ss_rate += exp( - pow(ss.q/R_/P->temperature, ss.stretch) ) *
//                     pow( von_mises/ss.sig0 , ss.n );

//  // Transient creep
//  AD::real F = 1;
//  for ( auto & tr : P->creep_md.tr )
//  {
//    AD::real etr_star = exp( tr.c * P->temperature ) * 
//                        pow( ( von_mises/tr.sig0 ), tr.m );

//    AD::real zeta = 0;
//    if ( abs(val(etr_star)) > 1e-20 ) zeta = 1 - P->creep_md.etr / etr_star;

//    double alpha = tr.alpha_w;
//    if ( zeta <= 0 ) alpha = 0;
//    F = F * exp( alpha * zeta * zeta );
//  }

//  AD::real etr_rate = ( F - 1 ) * creep_ss_rate;

//  P->creep_md.etr = P->creep_md.etr_n + val(etr_rate) * vpsolver.ts.dt;

//  /// NOTE: This must match the calculations in the update function (ViscoplasticIFC)
//  // Compute plastic strain rate
//  AD::Mat plastic_strain_rate = AD::Mat::Zero(3,3);

//  // If von_mises=0, the strain rate is zero
//  if ( von_mises )
//  for ( auto & ss : P->creep_md.ss )
//    plastic_strain_rate += 
//              3./2. * F * 
//              exp( - pow(ss.q/R_/P->temperature, ss.stretch) ) *
//              pow( von_mises/ss.sig0 , ss.n-1 ) *
//              deviatoric / ss.sig0;

//  // Update the plasteic strain
//  AD::Mat plastic_strain = vpsolver.ts.dt * plastic_strain_rate;
//  for (uint i=0; i<3; i++ )
//  for (uint j=0; j<3; j++ ) 
//    plastic_strain(i,j) += P->plastic_strain_n(i,j);

//  /// Set the results into the quadrature point
//  AD::dump( sigtot, P->sigtot );
//  AD::dump( sigeff, P->sigeff );
//  AD::dump( deviatoric, P->deviatoric );
//  AD::dump( plastic_strain_rate, P->plastic_strain_rate );
//  AD::dump( plastic_strain, P->plastic_strain);
//  AD::dump( grad_u, P->GRAD_U);

//  // Debug
////  const vector<Point> & xyz = fe->get_xyz();

//  P->von_mises = val(von_mises);
//  P->epskk = val(epskk);
//  P->F = val(F);
//  double g=0; 
//  for (uint i=0; i<3; i++ ) for (uint j=0; j<3; j++ ) g += pow( val(grad_u(i,j)), 2 );
//  P->grad_norm = g;

//  // 
//  // Subtract the plastic strain as a body force
//  //
//  // - ( \phi_j , Cijkl \varepsilon^p_kl ) 
//  for (uint B=0;  B<n_dofsv;  B++)
//  for (uint i=0; i<3; i++) 
//  for (uint j=0; j<3; j++) 
//  for (uint k=0; k<3; k++) 
//  for (uint l=0; l<3; l++) 
//    Fib(i,B) -= JxW[QP] * dphi[B][QP](j) * P->C_ijkl(i,j,k,l) * plastic_strain(k,l) ;

//  //   Hunting NaN
//  for (uint i=0; i<3; i++) 
//  for (uint B=0;  B<n_dofsv;  B++)
//    if ( P->plastic_strain.norm() > 1e30 ) 
//    {
//      ilog << "ELEM:" << elem->id();
//      ilog << "QP:" << QP;
//      ilog << "grad_u:" << Print(grad_u);
//      ilog << "sigeff:" << Print(sigeff);
//      ilog << "sigtot:" << Print(sigtot);
//      ilog << "deviatoric:" << Print(deviatoric);
//      ilog << "J2:" << J2;
//      ilog << "K:" << P->creep_md;
//      ilog << "T:" << P->temperature;
//      ilog << "VM:" << von_mises;
//      flog << "Found NaN in residual! " << *P;
//    }
// 
//  return ad_Fib;
}

/**
 *
 */
void VPMatEG::residual_and_jacobian_qp( EGFacePair & fp )
{
  // Lambda function to compatibilize stuff
  auto f = [this](const AD::Vec & x) { return this->residual_qp(x);  };

  Pp = fp.get_Pp(QP);
  Pn = fp.get_Pn(QP);

  AD::Vec F;
  ad.ad_Jijbm = AD::jacobian( f, wrt(ad.ad_Uib), at(ad.ad_Uib), F );

  // Update plastic_strain_k : the plastic strain after this newton iteration
  Pp->plastic_strain_k = Pp->plastic_strain;
  Pn->plastic_strain_k = Pn->plastic_strain;

  // Map from the AD variable to libmesh datastructures
  for (uint e=0; e<2; e++) 
  for (uint i=3; i<6; i++) 
  for (uint j=3; j<6; j++)
  for (uint B=0; B<Ndof(i);  B++)
  for (uint M=0; M<Ndof(j);  M++)
    Ke( ad.idx(e,i,B) , ad.idx(e,j,M) ) += val( ad.Jeijbm(e,i,j,B,M) );

  // Map from the AD variable to libmesh datastructures
  for (uint e=0; e<2; e++) 
  for (uint i=3; i<6; i++) 
  for (uint B=0;  B<Ndof(i);  B++)
    Re( ad.idx(e,i,B) ) += val( ad.Feib(e,i,B) );
}

/**
 *
 */
void VPMatEG::residual_and_jacobian ( const NumericVector<Number> & soln, 
                                      NumericVector<Number> * residual,
                                      SparseMatrix<Number> * jacobian )
{
  SCOPELOG(1);
  const DofMap & dof_map = system.get_dof_map();

  for ( EGFacePair & fp : gamma_I )
  {
    reinit( soln, fp );
    do { residual_and_jacobian_qp( fp ); } while ( next_qp() );
  }

  // Add to the residual 
  if ( residual ) {
    dof_map.constrain_element_vector ( Re, dof_indices_eg );
    residual->add_vector ( Re, dof_indices_eg );
  }

  // Add the the global matrix
  if ( jacobian ) {
    dof_map.constrain_element_matrix (Ke, dof_indices_eg);
    jacobian->add_matrix (Ke, dof_indices_eg);
  }
}

/**
 *
 */
using util::operator<<;
ostream& operator<<(ostream& os, const VPMatEG & m)
{
  os << "========= VPMatEG =======" << endl;
  os << "gamma_H (external boundary) : " << endl << m.gamma_H << endl;
  os << "gamma_I (internal skeleton) : " << endl << m.gamma_I << endl;
  os << "=========================" << endl;
  return os;
}

}} // ns
