#include "solver/viscoplastic/VPMatEG.h"

#include "config/MaterialConfig.h"
#include "solver/viscoplastic/VPSolver.h"

#include "libmesh/mesh.h"
#include "libmesh/dof_map.h"
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
fem_p(system), fem_n(system)
{ }

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
      if ( ! ( ( elem_n->level() <= elem_p->level() ) ||
             ( ( elem_n->active() ) &&
             ( elem_n->level() == elem_p->level() ) &&
             ( elem_p->id() < elem_n->id() ) ) ) ) continue; 
      
      EGFacePair egfp( elem_p->id(), side_p, elem_n->id() );

      fem_p.fe->reinit( elem_p, side_p );
      uint nqp = fem_p.qrule.n_points();
      dlog(1) << "==== NQP: " << nqp;

      egfp.Pq_p.resize(nqp);
      egfp.Pq_n.resize(nqp);

      gamma_I.emplace_back( egfp );
    }
    else
    {   // Register outer boundary skeleton
      EGFace egf( elem_p->id(), side_p );

      uint nqp = fem_p.qrule.n_points();
      egf.Pq.resize(nqp);

      gamma_H.emplace_back( egf );
    }
  }
}

}} // ns
