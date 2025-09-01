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
n_dofs(0),
n_dofsv(0),
n_uvars( vpsolver.is_eg() ? 6 : 3 ),
QP(0),
ad(n_uvars)
{ 
  // Init Qrule
  if ( ! system.has_variable( "UegX" ) )    flog << "Only EG systems should be here.";
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

  /* 3. Init QP counter and resolve autodiff solution vectors */
  QP = 0;
  // TODO: at some point the elements might have different number of dofs 
  ad.init( n_dofsv, n_dofs_eg, 2 ); 

  /* 4. Feed Uib with current solution */
  for ( uint e=0; e<2; e++ )
  for ( uint i=0; i<n_uvars; i++ )
  for ( uint B=0; B<Ndof(i); B++ )
    ad.Ueib(e,i,B) = soln( dof_indices[ad.idx(e,i,B)] );
}

/**
 *
 */
void VPMatEG::setup_dofs( EGFacePair & fp )
{
  MeshBase & mesh = system.get_mesh();
  Elem * elem_p = mesh.elem_ptr(fp.eid_p);
  Elem * elem_n = mesh.elem_ptr(fp.eid_n);

  const DofMap & dof_map = system.get_dof_map();
  vector<Elem*> ee = { elem_p, elem_n };

  fem_p.dof_indices.clear();
  fem_n.dof_indices.clear();

  /* Feed the dof indices in our datastructures */
  // 0..5 = "UX", "UY", "UZ", "UegX", "UegY", "UegZ"
  for ( uint vi=0; vi<6; vi++ )
  {
    vector<dof_id_type> dofi;

    dof_map.dof_indices ( elem_p, dofi, vi );
    fem_p.dof_indices.insert( fem_p.dof_indices.end(), dofi.begin(), dofi.end() );

    dof_map.dof_indices ( elem_n, dofi, vi );
    fem_n.dof_indices.insert( fem_n.dof_indices.end(), dofi.begin(), dofi.end() );
  }

  n_dofs  = fem_p.dof_indices.size();

  // Fill the unique dof_indices with all dofs
  dof_indices.clear();
  dof_indices.reserve( fem_p.dof_indices.size() + fem_n.dof_indices.size() );
  dof_indices.insert( dof_indices.end(), fem_p.dof_indices.begin() , fem_p.dof_indices.end() );
  dof_indices.insert( dof_indices.end(), fem_n.dof_indices.begin() , fem_n.dof_indices.end() );

  // 
  vector<dof_id_type> dof_indices_var;
  dof_map.dof_indices ( elem_p, dof_indices_var, 0 );
  n_dofsv = dof_indices_var.size();

  ASSERT( vpsolver.is_eg() , "Must be an EG solver to be here." );
  vector<dof_id_type> dof_indices_eg;
  dof_map.dof_indices (elem_p, dof_indices_eg, 3);
  n_dofs_eg = dof_indices_eg.size();

  /* Useful validation */
  dof_map.dof_indices (elem_n, dof_indices_eg, 3);
  ASSERT( n_dofs_eg == dof_indices_eg.size()  , "The elements of the interface have different number of EG DOFs (" << n_dofs_eg << " != " << dof_indices_eg.size() << ")");
  ASSERT( fp.Pq_p.size() == fp.Pq_n.size()    , "Size of the properties vector should be equal (" <<  fp.Pq_p.size() << " != " << fp.Pq_n.size() << ")" );
  ASSERT( fp.Pq_p.size() == qrule.n_points() , "Size of the properties vector should be equal to nqp (" << fp.Pq_p.size() << " != " << qrule.n_points() << ")" );


}

/**
 *
 */
void VPMatEG::residual_and_jacobian_qp()
{

}

/**
 *
 */
void VPMatEG::residual_and_jacobian ( const NumericVector<Number> & soln, 
                                      NumericVector<Number> * residual,
                                      SparseMatrix<Number> * jacobian )
{
  SCOPELOG(1);
  for ( EGFacePair & fp : gamma_I )
  {
    reinit( soln, fp );
    do { residual_and_jacobian_qp(); } while ( next_qp() );
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
