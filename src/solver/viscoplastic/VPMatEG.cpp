#include "solver/viscoplastic/VPMatEG.h"

#include "config/ModelConfig.h"
#include "config/MaterialConfig.h"
#include "solver/viscoplastic/VPSolver.h"
#include "util/OutputOperators.h"
#include "timeloop/Timestep.h"

#include "libmesh/mesh.h"
#include "libmesh/dof_map.h"
#include "libmesh/fe_map.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/dense_submatrix.h"
#include "libmesh/dense_subvector.h"
#include "libmesh/elem_side_builder.h"

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
ad(n_uvars),
elem_penalty(0)
{ 
  ASSERT( vpsolver.is_eg() , "Must be an EG solver to be here." );
  ASSERT( system.has_variable("UegX"), "Must be an EG solver to be here." );

  // Init Qrule based on a CG var (assuming it has higher order than EG)
  uint vid = system.variable_number( "UX" );
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

    const std::vector<Point> & xyz_p = fem_p.fe_cg->get_xyz();
    const std::vector<Point> & xyz_n = fem_n.fe_cg->get_xyz();
    for ( uint qp=0 ; qp<xyz_n.size() ; qp++ )
    {
      ASSERT( (xyz_p[qp]-xyz_n[qp]).norm() < 1e-8 , "These points should coincide.");
      const Point & pt = xyz_p[qp];
      fp.get_Pp(qp)->init_from_config( cfg_p, pt );
      fp.get_Pn(qp)->init_from_config( cfg_n, pt );
    }
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

      fem_p.fe_cg->reinit( elem_p, side_p );
      fem_p.fe_eg->reinit( elem_p, side_p );
      uint nqp = qrule.n_points();

      egfp.Pq_p.resize(nqp);
      egfp.Pq_n.resize(nqp);

      gamma_I.emplace_back( egfp );
    }
//    else
//    {   // Register outer boundary skeleton
//      EGFace egf( elem_p->id(), side_p );
//      uint nqp = qrule.n_points();
//      egf.Pq.resize(nqp);
//      gamma_D.emplace_back( egf );
//    }
  }

  init_properties();
}

/** **/
void VPMatEG::update_gammad()
{
  SCOPELOG(1);
  dlog(1) << "PROCESSING GAMMAD:";
  MeshBase & mesh = system.get_mesh();
  const BoundaryInfo & bi = mesh.get_boundary_info();

  // Boundary conditions : vpsolver.curr_bc.dirichlet => gamma_D
  const BC & curr_bc = vpsolver.curr_bc;
  for ( auto & dbc : curr_bc.dirichlet ) 
  {
    double val = dbc.val;
    boundary_id_type bid = dbc.bid;
    DirichletSetting diric;
    diric.vname = dbc.vname;
    diric.val = dbc.val;
    if ( ! system.has_variable( diric.vname ) ) flog << "System does not have variable '" << diric.vname << "'.";
    diric.vid_cg = system.variable_number( diric.vname );
    string vname_eg = diric.vname_eg();
    if ( ! system.has_variable( vname_eg ) ) flog << "System does not have variable '" << vname_eg << "'.";
    diric.vid_eg = system.variable_number( vname_eg );

    dlog(1) << "    DBC @ " << dbc.vname << " = " << dbc.val;

    for ( const auto & elem : mesh.active_local_element_ptr_range() )
    for ( auto side : elem->side_index_range() ) 
    if ( elem->neighbor_ptr(side) == nullptr )
    if ( bi.has_boundary_id( elem, side, bid ) )
    {
      EGFace egf( elem->id(), side );
      uint nqp = qrule.n_points();
      egf.Pq.resize(nqp);
      diric.egface_vec.push_back( egf );
    }
    gamma_D.push_back( diric );
  }
  using util::operator<<;
  dlog(1) << "GAMMAD:" << gamma_D;
}

/* Reinit structures for the face pair fp */
void VPMatEG::reinit( EGFacePair & fp )
{
  MeshBase & mesh = system.get_mesh();
  Elem * elem_p = mesh.elem_ptr(fp.eid_p);
  Elem * elem_n = mesh.elem_ptr(fp.eid_n);

  ElemSideBuilder side_builder;
  auto side_volume = side_builder(*elem_p, fp.side_p).volume();
//  uint elem_b_order = static_cast<unsigned int> (fe_elem_face->get_order());
  const Real h_elem = elem_p->volume()/side_volume; // * 1./pow(elem_b_order, 2.);
  elem_penalty = 1 / h_elem;

  /* 1. Init FEM structs */
  const std::vector<Point> & xyz_p = fem_p.fe_cg->get_xyz();
  fem_p.fe_cg->reinit( elem_p , fp.side_p );
  fem_p.fe_eg->reinit( elem_p , fp.side_p );

  vector<Point> pts;
  FEMap::inverse_map (3, elem_n, xyz_p, pts, 1E-10);
  fem_n.fe_cg->reinit(elem_n, & pts);
  fem_n.fe_eg->reinit(elem_n, & pts);
}

/* Reinit structures for the face pair fp */
void VPMatEG::reinit( const NumericVector<Number> & soln , EGFacePair & fp )
{
  /* 1. Basic reinit */
  reinit(fp);

  /* 2. Update dofmaps and dof counters. */
  setup_dofs( fp );

  uint nd = dof_indices.size();
  Ke.resize (nd, nd);
  Re.resize (nd);

  /* 3. Init QP counter and resolve autodiff solution vectors */
  QP = 0;


  // TODO: at some point the elements might have different number of dofs 
  //     : Only EG dofs
  ad.init( n_dofs_cg, n_dofs_eg, 2 ); 

  /* 4. Feed Uib with current solution - just EG part */
  for ( uint e=0; e<2; e++ )
  for ( uint i=0; i<6; i++ )
  for ( uint B=0; B<Ndof(i); B++ )
    ad.Ueib(e,i,B) = soln( dof_indices[ad.idx(e,i,B)] );

}

/**
 *
 */
void VPMatEG::setup_dofs( EGFacePair & fp )
{
  MeshBase & mesh = system.get_mesh();
  const DofMap & dof_map = system.get_dof_map();
  vector<dof_id_type> di;

  Elem * elem_p = mesh.elem_ptr(fp.eid_p);
  Elem * elem_n = mesh.elem_ptr(fp.eid_n);

  /* Feed the dof indices in our datastructures */
  fem_p.set_dofs( system, elem_p );
  fem_n.set_dofs( system, elem_n );

  // Flat structure with all dofs for both elements - for the dof assemblage
  dof_map.dof_indices ( elem_p, dof_indices );
  dof_map.dof_indices ( elem_n, di );
  dof_indices.insert( dof_indices.end(), di.begin() , di.end() );


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
  ad.ad_Fib.setZero(); 

  const Point & normal = fem_p.fe_cg->get_normals()[QP];
  const std::vector<Real> & JxW = fem_p.fe_cg->get_JxW();  

  // Compute U , jumps, averages
  const vector<vector<Real>> & phi_cg_p = fem_p.fe_cg->get_phi();
  const vector<vector<Real>> & phi_cg_n = fem_n.fe_cg->get_phi();
  const vector<vector<Real>> & phi_eg_p = fem_p.fe_eg->get_phi();
  const vector<vector<Real>> & phi_eg_n = fem_n.fe_eg->get_phi();
  AD::Vec u_p(3), u_n(3);
  // CG
  for (uint i=0; i<3; i++)
  for (uint M=0;  M<n_dofs_cg;  M++) 
  {
    u_p(i) += phi_cg_p[M][QP] * ad.Ueib(0,i,M);
    u_n(i) += phi_cg_n[M][QP] * ad.Ueib(1,i,M);
  }
  // EG
  for (uint i=0; i<3; i++)
  for (uint M=0;  M<n_dofs_eg;  M++) 
  {
    u_p(i) += phi_eg_p[M][QP] * ad.Ueib(0,i+3,M);
    u_n(i) += phi_eg_n[M][QP] * ad.Ueib(1,i+3,M);
  }
  AD::Vec jmp_u = u_p - u_n;
  dlog(1) << "jmp_u=" << jmp_u;

  // Compute grad_u
  const vector<vector<RealGradient>> & dphi_cg_p = fem_p.fe_cg->get_dphi();
  const vector<vector<RealGradient>> & dphi_cg_n = fem_n.fe_cg->get_dphi();
  const vector<vector<RealGradient>> & dphi_eg_p = fem_p.fe_eg->get_dphi();
  const vector<vector<RealGradient>> & dphi_eg_n = fem_n.fe_eg->get_dphi();
  AD::Mat grad_u_p(3,3), grad_u_n(3,3);
  // CG
  for (uint i=0; i<3; i++)
  for (uint j=0; j<3; j++) 
  for (uint M=0;  M<n_dofs_cg;  M++) 
  {
    grad_u_p(i, j) += dphi_cg_p[M][QP](j) * ad.Ueib(0,i,M);
    grad_u_n(i, j) += dphi_cg_n[M][QP](j) * ad.Ueib(1,i,M);
  }
  // EG
  for (uint i=0; i<3; i++)
  for (uint j=0; j<3; j++) 
  for (uint M=0;  M<n_dofs_eg;  M++) 
  {
    grad_u_p(i, j) += dphi_eg_p[M][QP](j) * ad.Ueib(0,i+3,M);
    grad_u_n(i, j) += dphi_eg_n[M][QP](j) * ad.Ueib(1,i+3,M);
  }
  dlog(1) << "grad_u_p:" << grad_u_p;
  dlog(1) << "grad_u_n:" << grad_u_n;

// Compute the stresses {sigma}n
//  sig = C_ijkl gradu_k,l
  AD::Mat sig_p(3,3), sig_n(3,3);
  for (uint i=0; i<3; i++) 
  for (uint j=0; j<3; j++) 
  for (uint k=0; k<3; k++) 
  for (uint l=0; l<3; l++) 
  {
    sig_p(i,j) += Pp->C_ijkl(i,j,k,l) * grad_u_p(k,l);
    sig_n(i,j) += Pn->C_ijkl(i,j,k,l) * grad_u_n(k,l);
  }
  for (uint i=0; i<3; i++) 
  for (uint j=0; j<3; j++) 
  {
    sig_p(i,j) += Pp->initial_stress(i,j);
    sig_n(i,j) += Pn->initial_stress(i,j);
  }
  for (uint k=0; k<3; k++ ) 
  {
    sig_p(k,k) -= Pp->alpha_d * ( Pp->temperature - Pp->initial_temperature );
    sig_p(k,k) -= Pp->biot * ( Pp->pressure - Pp->initial_pressure );
    sig_n(k,k) -= Pn->alpha_d * ( Pn->temperature - Pn->initial_temperature );
    sig_n(k,k) -= Pn->biot * ( Pn->pressure - Pn->initial_pressure );
  }

// Compute plasticity
  AD::Mat dev_p = sig_p, dev_n=sig_n;
  for (uint i=0; i<3; i++ )
  for (uint k=0; k<3; k++ ) 
  {
    dev_p(i,i) -= (1./3.) * sig_p(k,k);
    dev_n(i,i) -= (1./3.) * sig_n(k,k);
  }

  AD::real J2_p=0 , J2_n=0;
  for (uint i=0; i<3; i++ )
  for (uint j=0; j<3; j++ ) 
  {
    J2_p += (1./2.) * dev_p(i,j) * dev_p(i,j);
    J2_n += (1./2.) * dev_n(i,j) * dev_n(i,j);
  }

  AD::real von_mises_p = sqrt( 3 * J2_p );
  AD::real von_mises_n = sqrt( 3 * J2_n );
  double R_ = 8.3144;   // Universal gas constant [ J/mol/K ]

  AD::real creep_ss_rate_p = 0;
  for ( auto & ss : Pp->creep_md.ss )
    creep_ss_rate_p += exp( - pow(ss.q/R_/Pp->temperature, ss.stretch) ) *
                       pow( von_mises_p/ss.sig0 , ss.n );
  AD::real creep_ss_rate_n = 0;
  for ( auto & ss : Pn->creep_md.ss )
    creep_ss_rate_n += exp( - pow(ss.q/R_/Pn->temperature, ss.stretch) ) *
                       pow( von_mises_n/ss.sig0 , ss.n );

  // Transient creep (P)
  {
    AD::real F = 1;
    for ( auto & tr : Pp->creep_md.tr )
    {
      AD::real etr_star = exp( tr.c * Pp->temperature ) * 
        pow( ( von_mises_p/tr.sig0 ), tr.m );

      AD::real zeta = 0;
      if ( abs(val(etr_star)) > 1e-20 ) zeta = 1 - Pp->creep_md.etr / etr_star;

      double alpha = tr.alpha_w;
      if ( zeta <= 0 ) alpha = 0;
      F = F * exp( alpha * zeta * zeta );
    }
    AD::real etr_rate_p = ( F - 1 ) * creep_ss_rate_p;
    Pp->creep_md.etr = Pp->creep_md.etr_n + val(etr_rate_p) * vpsolver.ts.dt;

    AD::Mat plastic_strain_rate = AD::Mat::Zero(3,3);
    if ( von_mises_p )
    for ( auto & ss : Pp->creep_md.ss )
      plastic_strain_rate += 
                3./2. * F * 
                exp( - pow(ss.q/R_/Pp->temperature, ss.stretch) ) *
                pow( von_mises_p/ss.sig0 , ss.n-1 ) *
                dev_p / ss.sig0;

    AD::Mat plastic_strain = vpsolver.ts.dt * plastic_strain_rate;
    for (uint i=0; i<3; i++ )
    for (uint j=0; j<3; j++ ) 
      plastic_strain(i,j) += Pp->plastic_strain_n(i,j);

    // Update VPProps
    AD::dump( sig_p, Pp->sigtot );
    AD::dump( dev_p, Pp->deviatoric );
    AD::dump( plastic_strain_rate, Pn->plastic_strain_rate );
    AD::dump( plastic_strain, Pp->plastic_strain);
    AD::dump( grad_u_p, Pp->GRAD_U);
    Pp->von_mises = val(von_mises_p);

    // Update sig_p
    for (uint i=0; i<3; i++) 
    for (uint j=0; j<3; j++) 
    for (uint k=0; k<3; k++) 
    for (uint l=0; l<3; l++) 
      sig_p(i,j) -= Pp->C_ijkl(i,j,k,l) * plastic_strain(k,l);
  }

  // Transient creep (N)
  {
    AD::real F = 1;
    for ( auto & tr : Pn->creep_md.tr )
    {
      AD::real etr_star = exp( tr.c * Pn->temperature ) * 
        pow( ( von_mises_n/tr.sig0 ), tr.m );

      AD::real zeta = 0;
      if ( abs(val(etr_star)) > 1e-20 ) zeta = 1 - Pn->creep_md.etr / etr_star;

      double alpha = tr.alpha_w;
      if ( zeta <= 0 ) alpha = 0;
      F = F * exp( alpha * zeta * zeta );
    }
    AD::real etr_rate_n = ( F - 1 ) * creep_ss_rate_n;
    Pn->creep_md.etr = Pn->creep_md.etr_n + val(etr_rate_n) * vpsolver.ts.dt;

    AD::Mat plastic_strain_rate = AD::Mat::Zero(3,3);
    if ( von_mises_n )
    for ( auto & ss : Pp->creep_md.ss )
      plastic_strain_rate += 
                3./2. * F * 
                exp( - pow(ss.q/R_/Pn->temperature, ss.stretch) ) *
                pow( von_mises_n/ss.sig0 , ss.n-1 ) *
                dev_n / ss.sig0;

    AD::Mat plastic_strain = vpsolver.ts.dt * plastic_strain_rate;
    for (uint i=0; i<3; i++ )
    for (uint j=0; j<3; j++ ) 
      plastic_strain(i,j) += Pn->plastic_strain_n(i,j);

    // Update VPProps
    AD::dump( sig_n, Pn->sigtot );
    AD::dump( dev_n, Pn->deviatoric );
    AD::dump( plastic_strain_rate, Pn->plastic_strain_rate );
    AD::dump( plastic_strain, Pn->plastic_strain);
    AD::dump( grad_u_n, Pn->GRAD_U);
    Pn->von_mises = val(von_mises_n);

    // Update sig_n
    for (uint i=0; i<3; i++) 
    for (uint j=0; j<3; j++) 
    for (uint k=0; k<3; k++) 
    for (uint l=0; l<3; l++) 
      sig_n(i,j) -= Pn->C_ijkl(i,j,k,l) * plastic_strain(k,l);
  }

  // Calculate the jumps
  AD::Vec avg_sig_n = 0.5 * AD::dot( (sig_p+sig_n) , normal );
  AD::Vec jmp_sig_n =       AD::dot( (sig_p-sig_n) , normal );

  // Now we need to assemble Feib
  double beta = -1;
  double gamma = elem_penalty;

  // - ( [[w^B]] , {\sigma} n )_\GammaI
  for (uint B=0;  B<n_dofs_eg;  B++)
  for (uint i=0; i<3; i++) 
  for (uint j=0; j<3; j++) 
  {
    ad.Feib(0,i,B) -= JxW[QP] * phi_eg_p[B][QP] * avg_sig_n(i) ;
    ad.Feib(1,i,B) += JxW[QP] * phi_eg_n[B][QP] * avg_sig_n(i) ;
  }

  // + ( \gamma [[w^B]] , [[\sigma]] n )_\GammaI
  for (uint B=0;  B<n_dofs_eg;  B++)
  for (uint i=0; i<3; i++) 
  for (uint j=0; j<3; j++) 
  {
    ad.Feib(0,i,B) += JxW[QP] * gamma * phi_eg_p[B][QP] * jmp_sig_n(i) ;
    ad.Feib(1,i,B) -= JxW[QP] * gamma * phi_eg_n[B][QP] * jmp_sig_n(i) ;
  }

  // + ( \beta {\sigma_w}n , [[u]] )_\GammaI
  for (uint i=0; i<3; i++) 
  for (uint j=0; j<3; j++) 
  for (uint k=0; k<3; k++) 
  for (uint l=0; l<3; l++) 
  {
    // CG
    for (uint B=0;  B<n_dofs_cg;  B++)
    {
      ad.Feib(0,k,B) += 0.5 * JxW[QP] * beta * Pp->C_ijkl(i,j,k,l) * dphi_cg_p[B][QP](l) * normal(j) * jmp_u(i);
      ad.Feib(1,k,B) += 0.5 * JxW[QP] * beta * Pn->C_ijkl(i,j,k,l) * dphi_cg_n[B][QP](l) * normal(j) * jmp_u(i);
    }
    // EG
    // CG
    for (uint B=0;  B<n_dofs_eg;  B++)
    {
      ad.Feib(0,k,B) += 0.5 * JxW[QP] * beta * Pp->C_ijkl(i,j,k,l) * dphi_eg_p[B][QP](l) * normal(j);
      ad.Feib(1,k,B) += 0.5 * JxW[QP] * beta * Pn->C_ijkl(i,j,k,l) * dphi_eg_n[B][QP](l) * normal(j);
    }
  }

//  for (uint B=0;  B<n_dofsv;  B++)
//  for (uint i=0; i<3; i++) 
//  for (uint j=0; j<3; j++) 
//  for (uint k=0; k<3; k++) 
//  for (uint l=0; l<3; l++) 
//    Fib(i,B) += JxW[QP] *  dphi[B][QP](j) * P->C_ijkl(i,j,k,l) * grad_u(k,l) ;


  return ad.ad_Fib;

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
  for (uint i=0; i<6; i++) 
  for (uint j=0; j<6; j++)
  for (uint B=0; B<Ndof(i);  B++)
  for (uint M=0; M<Ndof(j);  M++)
    Ke( ad.idx(e,i,B) , ad.idx(e,j,M) ) += val( ad.Jeijbm(e,i,j,B,M) );

  // Map from the AD variable to libmesh datastructures
  for (uint e=0; e<2; e++) 
  for (uint i=0; i<6; i++) 
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
    dof_map.constrain_element_vector ( Re, dof_indices );
    residual->add_vector ( Re, dof_indices );
  }

  // Add the the global matrix
  if ( jacobian ) {
    dof_map.constrain_element_matrix (Ke, dof_indices);
    jacobian->add_matrix (Ke, dof_indices);
  }
}

/**
 *
 */
using util::operator<<;
ostream& operator<<(ostream& os, const VPMatEG & m)
{
  os << "========= VPMatEG =======" << endl;
  os << "gamma_D (external boundary) : " << endl << m.gamma_D << endl;
  os << "gamma_I (internal skeleton) : " << endl << m.gamma_I << endl;
  os << "=========================" << endl;
  return os;
}

ostream& operator<<(ostream& os, const DirichletSetting & m)
{
  os << endl << "Dirichlet Setting:";
  os << setw(20) << "vname:" << m.vname << "(" << m.vid_cg << "/" << m.vid_eg << ") = " << m.val;
  os << setw(20) << "egface_vec:" << endl;
  os << m.egface_vec;
  return os;
}

}} // ns
