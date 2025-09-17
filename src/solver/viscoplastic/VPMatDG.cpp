#include "solver/viscoplastic/VPMatDG.h"

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

VPMatDG::VPMatDG( TransientNonlinearImplicitSystem & sys_,
                  ViscoplasticSolver & vpsolver_ ) :
system(sys_),
vpsolver(vpsolver_),
fem_p(system),
fem_n(system),
qrule(2),
Pp(0), Pn(0),
n_uvars( 3 ),
n_dofs_u(0),
QP(0),
ad(n_uvars),
elem_penalty(0),
beta(-1)
{ 
  ASSERT( vpsolver.is_dg() , "Must be an DG solver to be here." );

  uint vid = system.variable_number( "UX" );
  DofMap & dof_map = system.get_dof_map();
  FEType fe_type = dof_map.variable_type(vid);
  dlog(1) << "QUAD_ORDER:" << fe_type.default_quadrature_order();
  qrule = QGauss( 2, fe_type.default_quadrature_order() );

  // Do not attach quadrature in fem_n -- it should be mapped by 
  // points during reinit (slave mapping)
  fem_p.attach_qrule( & qrule );
}

/**
 *
 */
void VPMatDG::init_properties()
{
  SCOPELOG(1);
  /* gamma_I */
  for ( DGFacePair & fp : gamma_I )
  {
    const MaterialConfig & cfg_p = vpsolver.get_material_config( fp.eid_p );
    const MaterialConfig & cfg_n = vpsolver.get_material_config( fp.eid_n );

    reinit_I( fp );

    const std::vector<Point> & xyz_p = fem_p.fe->get_xyz();
    const std::vector<Point> & xyz_n = fem_n.fe->get_xyz();
    for ( uint qp=0 ; qp<xyz_n.size() ; qp++ )
    {
      ASSERT( (xyz_p[qp]-xyz_n[qp]).norm() < 1e-8 , "These points should coincide.");
      const Point & pt = xyz_p[qp];
      fp.get_Pp(qp)->init_from_config( cfg_p, pt );
      fp.get_Pn(qp)->init_from_config( cfg_n, pt );
    }
  }
}

/** **/
void VPMatDG::update_variable( VARIABLE var , map<uint, double> val_by_sid )
{
  MeshBase & mesh = system.get_mesh();
  SCOPELOG(1);
  using util::operator<<;
  dlog(1) << "Updating variable " << v_to_str(var) << " from a val_by_sid map:" << val_by_sid;

  /** Gamma_I **/
  for ( DGFacePair & fp : gamma_I )
  {
    reinit_I( fp );
    do { 
      Elem * elem_p = mesh.elem_ptr(fp.eid_p);
      suint sid_p = elem_p->subdomain_id();

      Elem * elem_n = mesh.elem_ptr(fp.eid_n);
      suint sid_n = elem_n->subdomain_id();

      VPProps * Pp = fp.get_Pp(QP);
      VPProps * Pn = fp.get_Pn(QP);

      if ( ! val_by_sid.count( sid_p ) ) flog << "No initial value defined for element '" << fp.eid_p << "'.";
      if ( ! val_by_sid.count( sid_n ) ) flog << "No initial value defined for element '" << fp.eid_n << "'.";

      double val_p = val_by_sid[sid_p];
      double val_n = val_by_sid[sid_n];

      /** Update Pp and Pn**/
      if      ( var == TEMPERATURE )          { Pp->temperature = val_p;         Pn->temperature = val_n;          }
      else if ( var == INITIAL_TEMPERATURE )  { Pp->initial_temperature = val_p; Pn->initial_temperature = val_n;  }
      else if ( var == PRESSURE )             { Pp->pressure = val_p;            Pn->pressure = val_n;             }
      else if ( var == INITIAL_PRESSURE )     { Pp->initial_pressure = val_p;    Pn->initial_pressure = val_n;     }
      /** **/
    } while ( next_qp() );
  }

  /** Gamma_D **/
  for ( DGDirichlet & dgd : gamma_D )
  for ( DGFace & dgf : dgd.dgface_vec )
  {
    reinit_D( dgf );
    do {
      Elem * elem = mesh.elem_ptr(dgf.eid);
      suint sid = elem->subdomain_id();

      VPProps * P = dgf.get_P(QP);

      if ( ! val_by_sid.count( sid ) ) flog << "No initial value defined for element '" << elem << "'.";
      double val = val_by_sid[sid];

      /** Update Pp and Pn**/
      if      ( var == TEMPERATURE )          { P->temperature = val;          }
      else if ( var == INITIAL_TEMPERATURE )  { P->initial_temperature = val;  }
      else if ( var == PRESSURE )             { P->pressure = val;             }
      else if ( var == INITIAL_PRESSURE )     { P->initial_pressure = val;     }
      /** **/
    } while ( next_qp() );
  }

}

/**
 *
 */
void VPMatDG::init_properties_gammad()
{
  SCOPELOG(1);
  /* gamma_D */
  for ( DGDirichlet & dgd : gamma_D )
  {
    for ( DGFace & dgf : dgd.dgface_vec )
    {
      const MaterialConfig & cfg = vpsolver.get_material_config( dgf.eid );
      reinit_D( dgf );
      const std::vector<Point> & xyz = fem_p.fe->get_xyz();
      for ( uint qp=0 ; qp<xyz.size() ; qp++ )
        dgf.get_P(qp)->init_from_config( cfg, xyz[qp] );
    }
  }
}

/**
 *
 */
void VPMatDG::init()
{
  SCOPELOG(1);
  MeshBase & mesh = system.get_mesh();

  gamma_I.clear();

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
      
      DGFacePair dgfp( elem_p->id(), side_p, elem_n->id() );

      fem_p.fe->reinit( elem_p, side_p );
      uint nqp = qrule.n_points();

      dgfp.Pq_p.resize(nqp);
      dgfp.Pq_n.resize(nqp);

      gamma_I.emplace_back( dgfp );
    }
  }

  init_properties();
}

/** **/
void VPMatDG::update_gammad()
{
  SCOPELOG(1);

  gamma_D.clear();

  MeshBase & mesh = system.get_mesh();
  const BoundaryInfo & bi = mesh.get_boundary_info();

  // Boundary conditions : vpsolver.curr_bc.dirichlet => gamma_D
  const BC & curr_bc = vpsolver.curr_bc;
  for ( auto & dbc : curr_bc.dirichlet ) 
  {
    double val = dbc.val;
    boundary_id_type bid = dbc.bid;
    DGDirichlet diric;
    diric.vname = dbc.vname;
    diric.val = dbc.val;
    if ( ! system.has_variable( diric.vname ) ) flog << "System does not have variable '" << diric.vname << "'.";
    diric.vid = system.variable_number( diric.vname );

    for ( const auto & elem : mesh.active_local_element_ptr_range() )
    for ( auto side : elem->side_index_range() ) 
    if ( elem->neighbor_ptr(side) == nullptr )
    if ( bi.has_boundary_id( elem, side, bid ) )
    {
      DGFace dgf( elem->id(), side );
      reinit_D(dgf);
      uint nqp = qrule.n_points();
      dgf.Pq.resize(nqp);
      diric.dgface_vec.push_back( dgf );
    }
    gamma_D.push_back( diric );
  }

  /** TODO : In this architecture, we are gonna miss the plasticity at the interface! **/
  init_properties_gammad();

  using util::operator<<;
  dlog(1) << "GAMMAD:" << gamma_D;
}

/* Reinit structures for the face pair fp */
void VPMatDG::reinit_I( DGFacePair & fp )
{
  MeshBase & mesh = system.get_mesh();
  Elem * elem_p = mesh.elem_ptr(fp.eid_p);
  Elem * elem_n = mesh.elem_ptr(fp.eid_n);

  ElemSideBuilder side_builder;
  auto side_volume = side_builder(*elem_p, fp.side_p).volume();
//  uint elem_b_order = static_cast<unsigned int> (fe_elem_face->get_order());
  const Real h_elem = elem_p->volume()/side_volume; // * 1./pow(elem_b_order, 2.);
//  dlog(1) << "h_elem:" << h_elem;
  elem_penalty = 1e12 / h_elem;

  /* 1. Init FEM structs */
  const std::vector<Point> & xyz_p = fem_p.fe->get_xyz();
  fem_p.fe->reinit( elem_p , fp.side_p );

  vector<Point> pts;
  FEMap::inverse_map (3, elem_n, xyz_p, pts, 1E-10);
  fem_n.fe->reinit(elem_n, & pts);

  QP = 0;
}

/* Reinit structures for the face pair fp */
void VPMatDG::reinit_I( const NumericVector<Number> & soln , DGFacePair & fp )
{
  /* 1. Basic reinit */
  reinit_I(fp);

  /* 2. Update dofmaps and dof counters. */
  setup_dofs( fp );

  uint nd = dof_indices.size();
  Ke.resize (nd, nd);
  Re.resize (nd);

  ad.init( n_dofs_u, 2 ); 

  /* 4. Feed Uib with current solution - just DG part */
  for ( uint e=0; e<2; e++ )
  for ( uint i=0; i<3; i++ )
  for ( uint B=0; B<n_dofs_u; B++ )
    ad.Ueib(e,i,B) = soln( dof_indices[ad.idx(e,i,B)] );
}

/* Reinit structures for the face pair fp */
void VPMatDG::reinit_D( DGFace & dgf )
{
  SCOPELOG(1);
  MeshBase & mesh = system.get_mesh();
  Elem * elem = mesh.elem_ptr(dgf.eid);

  ElemSideBuilder side_builder;
  auto side_volume = side_builder(*elem, dgf.side).volume();
//  uint elem_b_order = static_cast<unsigned int> (fe_elem_face->get_order());
  const Real h_elem = elem->volume()/side_volume; // * 1./pow(elem_b_order, 2.);
//  dlog(1) << "h_elem:" << h_elem;
  elem_penalty = 1e12 / h_elem;

  /* Init FEM structs */
  fem_p.fe->reinit( elem , dgf.side );

  QP = 0;
}
/* Reinit structures for the face pair fp */
void VPMatDG::reinit_D( const NumericVector<Number> & soln , DGFace & dgf )
{
  /* 1. Basic reinit */
  reinit_D(dgf);

  /* 2. Update dofmaps and dof counters. */
  setup_dofs( dgf );

  uint nd = dof_indices.size();
  Ke.resize (nd, nd);
  Re.resize (nd);

  ad.init( n_dofs_u, 1 ); 

  /* 4. Feed Uib with current solution */
  for ( uint i=0; i<3; i++ )
  for ( uint B=0; B<n_dofs_u; B++ )
    ad.Uib(i,B) = soln( dof_indices[ad.idx(i,B)] );
}

/**
 *
 */
void VPMatDG::setup_dofs( DGFace & dgf )
{
  MeshBase & mesh = system.get_mesh();
  const DofMap & dof_map = system.get_dof_map();
  vector<dof_id_type> di;

  Elem * elem = mesh.elem_ptr(dgf.eid);

  /* Feed the dof indices in our datastructures */
  fem_p.set_dofs( system, elem );

  // Flat structure with all dofs for both elements - for the dof assemblage
  dof_map.dof_indices ( elem, dof_indices );

  dof_map.dof_indices ( elem, di, 0 );
  n_dofs_u = di.size();
}

/**
 *
 */
void VPMatDG::setup_dofs( DGFacePair & fp )
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
  n_dofs_u = di.size();

  /* Useful validation */
  ASSERT( fp.Pq_p.size() == fp.Pq_n.size()    , "Size of the properties vector should be equal (" <<  fp.Pq_p.size() << " != " << fp.Pq_n.size() << ")" );
  ASSERT( fp.Pq_p.size() == qrule.n_points() , "Size of the properties vector should be equal to nqp (" << fp.Pq_p.size() << " != " << qrule.n_points() << ")" );
}

/**
 * 
 * This is the actual calculation of the residual using the AutoDiff types.
 *
 */
AD::Vec VPMatDG::residual_gamma_I_qp( const AD::Vec & /* ad_Uib */ )
{
  ASSERT( Pp->lame_mu || Pn->lame_mu, "Null lame_mu?" );
  ad.ad_Fib.setZero(); 

  const Point & normal = fem_p.fe->get_normals()[QP];
  const std::vector<Real> & JxW = fem_p.fe->get_JxW();  

  // Compute U , jumps, averages
  const vector<vector<Real>> & phi_p = fem_p.fe->get_phi();
  const vector<vector<Real>> & phi_n = fem_n.fe->get_phi();
  AD::Vec u_p(3), u_n(3);
  // CG
  for (uint i=0; i<3; i++)
  for (uint M=0;  M<n_dofs_u;  M++) 
  {
    u_p(i) += phi_p[M][QP] * ad.Ueib(0,i,M);
    u_n(i) += phi_n[M][QP] * ad.Ueib(1,i,M);
  }

  AD::Vec jmp_u = u_p - u_n;

  // Compute grad_u
  const vector<vector<RealGradient>> & dphi_p = fem_p.fe->get_dphi();
  const vector<vector<RealGradient>> & dphi_n = fem_n.fe->get_dphi();
  AD::Mat grad_u_p(3,3), grad_u_n(3,3);
  // CG
  for (uint i=0; i<3; i++)
  for (uint j=0; j<3; j++) 
  for (uint M=0;  M<n_dofs_u;  M++) 
  {
    grad_u_p(i, j) += dphi_p[M][QP](j) * ad.Ueib(0,i,M);
    grad_u_n(i, j) += dphi_n[M][QP](j) * ad.Ueib(1,i,M);
  }

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
//  for (uint k=0; k<3; k++ ) 
//  {
//    sig_p(k,k) -= Pp->alpha_d * ( Pp->temperature - Pp->initial_temperature );
//    sig_p(k,k) -= Pp->biot * ( Pp->pressure - Pp->initial_pressure );
//    sig_n(k,k) -= Pn->alpha_d * ( Pn->temperature - Pn->initial_temperature );
//    sig_n(k,k) -= Pn->biot * ( Pn->pressure - Pn->initial_pressure );
//  }

//// Compute plasticity
//  AD::Mat dev_p = sig_p, dev_n=sig_n;
//  for (uint i=0; i<3; i++ )
//  for (uint k=0; k<3; k++ ) 
//  {
//    dev_p(i,i) -= (1./3.) * sig_p(k,k);
//    dev_n(i,i) -= (1./3.) * sig_n(k,k);
//  }

//  AD::real J2_p=0 , J2_n=0;
//  for (uint i=0; i<3; i++ )
//  for (uint j=0; j<3; j++ ) 
//  {
//    J2_p += (1./2.) * dev_p(i,j) * dev_p(i,j);
//    J2_n += (1./2.) * dev_n(i,j) * dev_n(i,j);
//  }

//  AD::real von_mises_p = sqrt( 3 * J2_p );
//  AD::real von_mises_n = sqrt( 3 * J2_n );
//  double R_ = 8.3144;   // Universal gas constant [ J/mol/K ]

//  AD::real creep_ss_rate_p = 0;
//  for ( auto & ss : Pp->creep_md.ss )
//    creep_ss_rate_p += exp( - pow(ss.q/R_/Pp->temperature, ss.stretch) ) *
//                       pow( von_mises_p/ss.sig0 , ss.n );
//  AD::real creep_ss_rate_n = 0;
//  for ( auto & ss : Pn->creep_md.ss )
//    creep_ss_rate_n += exp( - pow(ss.q/R_/Pn->temperature, ss.stretch) ) *
//                       pow( von_mises_n/ss.sig0 , ss.n );

  // Transient creep (P)
//  {
//    AD::real F = 1;
//    for ( auto & tr : Pp->creep_md.tr )
//    {
//      AD::real etr_star = exp( tr.c * Pp->temperature ) * 
//        pow( ( von_mises_p/tr.sig0 ), tr.m );

//      AD::real zeta = 0;
//      if ( abs(val(etr_star)) > 1e-20 ) zeta = 1 - Pp->creep_md.etr / etr_star;

//      double alpha = tr.alpha_w;
//      if ( zeta <= 0 ) alpha = 0;
//      F = F * exp( alpha * zeta * zeta );
//    }
//    AD::real etr_rate_p = ( F - 1 ) * creep_ss_rate_p;
//    Pp->creep_md.etr = Pp->creep_md.etr_n + val(etr_rate_p) * vpsolver.ts.dt;

//    AD::Mat plastic_strain_rate = AD::Mat::Zero(3,3);
//    if ( von_mises_p )
//    for ( auto & ss : Pp->creep_md.ss )
//      plastic_strain_rate += 
//                3./2. * F * 
//                exp( - pow(ss.q/R_/Pp->temperature, ss.stretch) ) *
//                pow( von_mises_p/ss.sig0 , ss.n-1 ) *
//                dev_p / ss.sig0;

//    AD::Mat plastic_strain = vpsolver.ts.dt * plastic_strain_rate;
//    for (uint i=0; i<3; i++ )
//    for (uint j=0; j<3; j++ ) 
//      plastic_strain(i,j) += Pp->plastic_strain_n(i,j);

//    // Update VPProps
//    AD::dump( sig_p, Pp->sigtot );
//    AD::dump( dev_p, Pp->deviatoric );
//    AD::dump( plastic_strain_rate, Pn->plastic_strain_rate );
//    AD::dump( plastic_strain, Pp->plastic_strain);
//    AD::dump( grad_u_p, Pp->GRAD_U);
//    Pp->von_mises = val(von_mises_p);

//    // Update sig_p
//    for (uint i=0; i<3; i++) 
//    for (uint j=0; j<3; j++) 
//    for (uint k=0; k<3; k++) 
//    for (uint l=0; l<3; l++) 
//      sig_p(i,j) -= Pp->C_ijkl(i,j,k,l) * plastic_strain(k,l);
//  }

  // Transient creep (N)
//  {
//    AD::real F = 1;
//    for ( auto & tr : Pn->creep_md.tr )
//    {
//      AD::real etr_star = exp( tr.c * Pn->temperature ) * 
//        pow( ( von_mises_n/tr.sig0 ), tr.m );

//      AD::real zeta = 0;
//      if ( abs(val(etr_star)) > 1e-20 ) zeta = 1 - Pn->creep_md.etr / etr_star;

//      double alpha = tr.alpha_w;
//      if ( zeta <= 0 ) alpha = 0;
//      F = F * exp( alpha * zeta * zeta );
//    }
//    AD::real etr_rate_n = ( F - 1 ) * creep_ss_rate_n;
//    Pn->creep_md.etr = Pn->creep_md.etr_n + val(etr_rate_n) * vpsolver.ts.dt;

//    AD::Mat plastic_strain_rate = AD::Mat::Zero(3,3);
//    if ( von_mises_n )
//    for ( auto & ss : Pp->creep_md.ss )
//      plastic_strain_rate += 
//                3./2. * F * 
//                exp( - pow(ss.q/R_/Pn->temperature, ss.stretch) ) *
//                pow( von_mises_n/ss.sig0 , ss.n-1 ) *
//                dev_n / ss.sig0;

//    AD::Mat plastic_strain = vpsolver.ts.dt * plastic_strain_rate;
//    for (uint i=0; i<3; i++ )
//    for (uint j=0; j<3; j++ ) 
//      plastic_strain(i,j) += Pn->plastic_strain_n(i,j);

//    // Update VPProps
//    AD::dump( sig_n, Pn->sigtot );
//    AD::dump( dev_n, Pn->deviatoric );
//    AD::dump( plastic_strain_rate, Pn->plastic_strain_rate );
//    AD::dump( plastic_strain, Pn->plastic_strain);
//    AD::dump( grad_u_n, Pn->GRAD_U);
//    Pn->von_mises = val(von_mises_n);

//    // Update sig_n
//    for (uint i=0; i<3; i++) 
//    for (uint j=0; j<3; j++) 
//    for (uint k=0; k<3; k++) 
//    for (uint l=0; l<3; l++) 
//      sig_n(i,j) -= Pn->C_ijkl(i,j,k,l) * plastic_strain(k,l);
//  }

  // Calculate the jumps
  // The normal is outward to the "p" face, so we need the negative sign for sig_n
//  AD::Vec avg_sig_n = 0.5 * AD::dot( (sig_p + sig_n) , normal );

  // Now we need to assemble Feib
  double gamma = 1e15; //elem_penalty;

  // - ( [[w^B]] , {\sigma n} )_\GammaI
  for (uint B=0;  B<n_dofs_u;  B++)
  for (uint i=0; i<3; i++) 
  for (uint j=0; j<3; j++) 
  {
    ad.Feib(0,i,B) -= 0.5 * JxW[QP] * phi_p[B][QP] * ( sig_p(i,j) + sig_n(i,j) ) * normal(j) ;
    ad.Feib(1,i,B) += 0.5 * JxW[QP] * phi_n[B][QP] * ( sig_p(i,j) + sig_n(i,j) ) * normal(j) ;
  }

  // + ( \gamma [[w^B]] , [[u]] )_\GammaI
  for (uint B=0;  B<n_dofs_u;  B++)
  for (uint i=0; i<3; i++) 
  {
    ad.Feib(0,i,B) += JxW[QP] * gamma * phi_p[B][QP] * jmp_u(i) ;
    ad.Feib(1,i,B) -= JxW[QP] * gamma * phi_n[B][QP] * jmp_u(i) ;
  }

  // + ( \beta {\sigma_w} n , [[u]] )_\GammaI
  for (uint i=0; i<3; i++) 
  for (uint j=0; j<3; j++) 
  for (uint k=0; k<3; k++) 
  for (uint l=0; l<3; l++) 
  {
    for (uint B=0;  B<n_dofs_u;  B++)
    {
      // normal(j) is the outward normal to the "p" face ; The "n" normal requires the inverse sign.
      ad.Feib(0,k,B) += 0.5 * JxW[QP] * beta * Pp->C_ijkl(i,j,k,l) * dphi_p[B][QP](l) * normal(j) * jmp_u(i); // jmp_u=0 for CG?
      ad.Feib(1,k,B) += 0.5 * JxW[QP] * beta * Pn->C_ijkl(i,j,k,l) * dphi_n[B][QP](l) * normal(j) * jmp_u(i);
    }
  }

//  // Do we need this?
//  for (uint i=0; i<3; i++) 
//  for (uint j=0; j<3; j++) 
//  for (uint B=0;  B<n_dofs_u;  B++)
//  {
//    ad.Feib(0,i,B) -= 0.5 * JxW[QP] * beta * Pp->alpha_d * ( Pp->temperature - Pp->initial_temperature ) * normal(j) * jmp_u(i); 
//    ad.Feib(0,i,B) -= 0.5 * JxW[QP] * beta * Pp->biot    * ( Pp->pressure    - Pp->initial_pressure    ) * normal(j) * jmp_u(i);
//    ad.Feib(1,i,B) -= 0.5 * JxW[QP] * beta * Pn->alpha_d * ( Pn->temperature - Pn->initial_temperature ) * normal(j) * jmp_u(i); 
//    ad.Feib(1,i,B) -= 0.5 * JxW[QP] * beta * Pn->biot    * ( Pn->pressure    - Pn->initial_pressure    ) * normal(j) * jmp_u(i);
//  }

  return ad.ad_Fib;
}

/**
 *
 */
void VPMatDG::residual_and_jacobian_gammaI_qp( DGFacePair & fp )
{
  // Lambda function to compatibilize stuff
  auto f = [this](const AD::Vec & x) { return this->residual_gamma_I_qp(x);  };

  Pp = fp.get_Pp(QP);
  Pn = fp.get_Pn(QP);

  if ( Pp->temperature < 390 ) dlog(1) << "Pp->temperature: " << Pp->temperature;
  if ( Pn->temperature < 390 ) dlog(1) << "Pn->temperature: " << Pn->temperature;

  AD::Vec F;
  ad.ad_Jijbm = AD::jacobian( f, wrt(ad.ad_Uib), at(ad.ad_Uib), F );

  // Update plastic_strain_k : the plastic strain after this newton iteration
  Pp->plastic_strain_k = Pp->plastic_strain;
  Pn->plastic_strain_k = Pn->plastic_strain;

  // Map from the AD variable to libmesh datastructures
  for (uint e=0; e<2; e++) 
  for (uint n=0; n<2; n++) 
  for (uint i=0; i<3; i++) 
  for (uint j=0; j<3; j++)
  for (uint B=0; B<n_dofs_u;  B++)
  for (uint M=0; M<n_dofs_u;  M++)
    Ke( ad.idx(e,i,B) , ad.idx(n,j,M) ) += val( ad.Jenijbm(e,n,i,j,B,M) );

  // Map from the AD variable to libmesh datastructures
  for (uint e=0; e<2; e++) 
  for (uint i=0; i<3; i++) 
  for (uint B=0;  B<n_dofs_u;  B++)
    Re( ad.idx(e,i,B) ) += val( ad.Feib(e,i,B) );
}

/**
 *
 */
void VPMatDG::residual_and_jacobian ( const NumericVector<Number> & soln, 
                                      NumericVector<Number> * residual,
                                      SparseMatrix<Number> * jacobian )
{
  SCOPELOG(1);
  const DofMap & dof_map = system.get_dof_map();

//  dlog(1) << "GAMMA_I" << gamma_I;

  /* 
   * Gamma_I 
   */
  for ( DGFacePair & fp : gamma_I )
  {
    reinit_I( soln, fp );
    do { residual_and_jacobian_gammaI_qp( fp ); } while ( next_qp() );

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

  /*
   * Gamma_D
   */
  for ( DGDirichlet & dgd : gamma_D )
  {
    /* The dirichlet value we are setting */
    u_hat = dgd.u_hat();
    for ( DGFace & dgf : dgd.dgface_vec )
    {
      reinit_D( soln, dgf );

      do { residual_and_jacobian_gammaD_qp( dgf ); } while ( next_qp() );

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
  }
}

/**
 * 
 * This is the actual calculation of the residual using the AutoDiff types.
 *
 */
AD::Vec VPMatDG::residual_gammaD_qp( const AD::Vec & /* ad_Uib */ )
{
  ASSERT( Pp->lame_mu , "Null lame_mu?" );
  ad.ad_Fib.setZero(); 

  const Point & normal = fem_p.fe->get_normals()[QP];
  const std::vector<Real> & JxW = fem_p.fe->get_JxW();  

  // Compute U , jumps, averages
  const vector<vector<Real>> & phi = fem_p.fe->get_phi();
  const vector<vector<RealGradient>> & dphi = fem_p.fe->get_dphi();

  AD::Vec u(3);

  // CG
  for (uint i=0; i<3; i++)
  for (uint M=0;  M<n_dofs_u;  M++) 
    u(i) += phi[M][QP] * ad.Ueib(0,i,M);

  AD::Mat grad_u(3,3);
  // CG
  for (uint i=0; i<3; i++)
  for (uint j=0; j<3; j++) 
  for (uint M=0;  M<n_dofs_u;  M++) 
    grad_u(i, j) += dphi[M][QP](j) * ad.Ueib(0,i,M);

  // Now we need to assemble Feib
  double gamma = 1e15; //elem_penalty;

  // + ( \gamma w^B , u - u_hat n )_\GammaD
  for (uint i=0; i<3; i++) 
  if (u_hat[i]) 
  for (uint B=0;  B<n_dofs_u;  B++)      // CG
    ad.Fib(i,B) += JxW[QP] * gamma * phi[B][QP] * ( u(i) - *(u_hat[i]) ) ;

  // + ( \beta {\sigma_w} n , u - u_hat )_\GammaD
  for (uint i=0; i<3; i++) 
  if (u_hat[i]) 
  for (uint j=0; j<3; j++) 
  for (uint k=0; k<3; k++) 
  if (u_hat[k]) 
  for (uint l=0; l<3; l++) 
  for (uint B=0;  B<n_dofs_u;  B++)
    ad.Fib(k,B) += JxW[QP] * beta * Pp->C_ijkl(i,j,k,l) * dphi[B][QP](l) * normal(j) * ( u(i) - *(u_hat[i]) );

//  for (uint i=0; i<3; i++) 
//  if (u_hat[i]) 
//  for (uint j=0; j<3; j++) 
//  for (uint B=0;  B<n_dofs_u;  B++)
//  {
//    ad.Fib(i,B) -= JxW[QP] * beta * Pp->alpha_d * ( Pp->temperature - Pp->initial_temperature ) * normal(j) * ( u(i) - *(u_hat[i]) );
//    ad.Fib(i,B) -= JxW[QP] * beta * Pp->biot *    ( Pp->pressure    - Pp->initial_pressure )    * normal(j) * ( u(i) - *(u_hat[i]) );
//  }

  // - ( v  , {\sigma_u} n )_\GammaD
  for (uint i=0; i<3; i++) 
  if  ( u_hat[i] ) 
  for (uint j=0; j<3; j++) 
  for (uint k=0; k<3; k++) 
  if (u_hat[k]) 
  for (uint l=0; l<3; l++) 
  for (uint B=0;  B<n_dofs_u;  B++)
    ad.Fib(i,B) -= JxW[QP] *  Pp->C_ijkl(i,j,k,l) * phi[B][QP] * normal(j) * grad_u(k,l) ;

//  for (uint i=0; i<3; i++) 
//  if  ( u_hat[i] ) 
//  for (uint j=0; j<3; j++) 
//  for (uint B=0;  B<n_dofs_u;  B++)
//  {
//    ad.Fib(i,B) -= JxW[QP] * phi[B][QP] * normal(j) * Pp->alpha_d * ( Pp->temperature - Pp->initial_temperature ) ;
//    ad.Fib(i,B) -= JxW[QP] * phi[B][QP] * normal(j) * Pp->biot    * ( Pp->pressure    - Pp->initial_pressure ) ;
//  }

  return ad.ad_Fib;
}

/**
 *
 */
void VPMatDG::residual_and_jacobian_gammaD_qp( DGFace & dgf )
{
  // Lambda function to compatibilize stuff
  auto f = [this](const AD::Vec & x) { return this->residual_gammaD_qp(x);  };

  dlog(1) << "residual_and_jacobian_gammaD_qp : QP:" << QP ;
  Pp = dgf.get_P(QP);

  AD::Vec F;
  ad.ad_Jijbm = AD::jacobian( f, wrt(ad.ad_Uib), at(ad.ad_Uib), F );

  // Map from the AD variable to libmesh datastructures
  for (uint i=0; i<3; i++) 
  for (uint j=0; j<3; j++)
  for (uint B=0; B<n_dofs_u;  B++)
  for (uint M=0; M<n_dofs_u;  M++)
    Ke( ad.idx(i,B) , ad.idx(j,M) ) += val( ad.Jijbm(i,j,B,M) );

  // Map from the AD variable to libmesh datastructures
  for (uint i=0; i<3; i++) 
  for (uint B=0;  B<n_dofs_u;  B++)
    Re( ad.idx(i,B) ) += val( ad.Fib(i,B) );

//  using util::operator<<;
//  dlog(1) << "Re: " << Re;
//  dlog(1) << "Ke: " << Ke;
}

 /**
 *
 */
using util::operator<<;
ostream& operator<<(ostream& os, const VPMatDG & m)
{
  os << "========= VPMatDG =======" << endl;
  os << "gamma_D (external boundary) : " << endl << m.gamma_D << endl;
  os << "gamma_I (internal skeleton) : " << endl << m.gamma_I << endl;
  os << "=========================" << endl;
  return os;
}

ostream& operator<<(ostream& os, const DGDirichlet & m)
{
  os << endl << "Dirichlet Setting:";
  os << setw(20) << "vname:" << m.vname << "(" << m.vid << ") = " << m.val;
  os << setw(20) << "dgface_vec:" << endl;
  os << m.dgface_vec;
  return os;
}

}} // ns
