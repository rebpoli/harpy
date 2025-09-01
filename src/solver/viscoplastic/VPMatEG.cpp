#include "solver/viscoplastic/VPMatEG.h"

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
n_dofsv(0),
n_uvars( vpsolver.is_eg() ? 6 : 3 ),
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
}

/* Reinit structures for the face pair fp */
void VPMatEG::reinit( const NumericVector<Number> & soln , EGFacePair & fp )
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
  {
    dlog(1) << "n_dofs_eg:" << n_dofs_eg << " e:" << e << " i:" << i << " B" << B  << " ad.idx:" << ad.idx(e,i,B);
    ad.Ueib(e,i,B) = soln( dof_indices_eg[ad.idx(e,i,B)] );
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

  // Init dofs counters
  const DofMap & dof_map = system.get_dof_map();
  vector<dof_id_type> di;
  dof_map.dof_indices ( elem_p, di, 0 );
  n_dofsv = di.size();
  dof_map.dof_indices ( elem_p, di, 3 );
  n_dofs_eg = di.size();

  /* Useful validation */
  ASSERT( fp.Pq_p.size() == fp.Pq_n.size()    , "Size of the properties vector should be equal (" <<  fp.Pq_p.size() << " != " << fp.Pq_n.size() << ")" );
  ASSERT( fp.Pq_p.size() == qrule.n_points() , "Size of the properties vector should be equal to nqp (" << fp.Pq_p.size() << " != " << qrule.n_points() << ")" );
}

/**
 *
 */
void VPMatEG::residual_and_jacobian_qp()
{
  /* */

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
    do { residual_and_jacobian_qp(); } while ( next_qp() );
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
