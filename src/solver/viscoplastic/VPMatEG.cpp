#include "solver/viscoplastic/VPMatEG.h"

#include "config/MaterialConfig.h"
#include "solver/viscoplastic/VPSolver.h"
#include "util/OutputOperators.h"

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
fem_p(system),
fem_n(system),
qrule(2)
{ 
  // Init Qrule
  if ( ! system.has_variable( "UegX" ) )    flog << "Only EG systems should be here.";
  uint vid = system.variable_number( "UegX" );
  DofMap & dof_map = system.get_dof_map();
  FEType fe_type = dof_map.variable_type(vid);
  qrule = QGauss( 2, fe_type.default_quadrature_order() );
  // Do not attach quadrature in fem_n -- it should be mapped by points (slave mapping)
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
      // Reinit
//      reinit( soln, fp )
//      do { residual_and_jacobian_qp(); } while ( next_qp() );
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
