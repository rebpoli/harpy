
#include "base/Global.h"
#include "solver/ElemProjection.h"
#include "util/Stopwatch.h"

#include "libmesh/system.h"
#include "libmesh/dof_map.h"
#include "libmesh/elem.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/explicit_system.h"
#include "libmesh/numeric_vector.h"

using namespace libMesh;

/**
 *
 *
 */
ElemProjection::ElemProjection( ExplicitSystem & sys_, unique_ptr<FEBase> & fe_, QGauss * qrule_ ) : 
    sys(sys_), fe(fe_), qrule(qrule_), elem(0) 
{
  // Importante reservar as variaveis que vamos usar na avaliacao, para
  // elas serem utilizadas no reinit, mesmo que o reinit seja feito externamente a essa classe
  fe->get_phi();
  fe->get_JxW();
  fe->get_xyz();
}

/**
 *
 *
 *
 */
void ElemProjection::reinit( Elem * e, bool fe_reinit, int side ) 
{
  elem = e;
  if ( fe_reinit ) {
    if ( side < 0 ) fe->reinit(elem);
    else fe->reinit(elem, side);
  }
}

/**
 *
 *
 *
 */
void ElemProjection::eval( string vname, const vector<double> & vals_qp ) {
  dlog(5) << "Projetando variavel '"<<vname<<"' ("<< sys.variable_number(vname)<<")...";
  sys.print_info();
  eval( sys.variable_number(vname), vals_qp );
}

/**
 *
 *
 */
void ElemProjection::eval(uint var, const vector<double> & vals_qp )
{
  if ( ! elem ) flog << "Elemento nulo! Como avaliar?";

  const std::vector<std::vector<Real>> & phi = fe->get_phi();
  const std::vector<Real> & jxw = fe->get_JxW();

  // Matriz e vetor do sistema local do elemento: M . x = F
  uint n_dofs = phi.size();
  DenseMatrix<Real> M(n_dofs, n_dofs);
  DenseVector<Number> F(n_dofs);

  // Monta matriz
//  dlog(1) << "Montando matrizes ... (nqp="<<qrule->n_points()<<") ndofs:" << n_dofs;
  for (uint qp=0; qp<qrule->n_points(); qp++) {
//    dlog(1) << "xyz["<<qp<<"] : " << xyz[qp];
    for(uint i=0; i<n_dofs; i++) 
    {
      F(i) += jxw[qp] * vals_qp[qp] * phi[i][qp];

      for(uint j=0; j<n_dofs; j++)
        M(i,j) += jxw[qp] * phi[i][qp] * phi[j][qp];
    }

    for(uint i=0; i<n_dofs; i++) 
      if ( std::abs(M(i,i)) < 1e-10 ) M(i,i) = 1;
  }

//  dlog(1) << "Solving: F=" << F;
//  dlog(1) << "M=" << M;
  DenseVector<Number> projected_data;
//  M.cholesky_solve(F, projected_data);
  M.lu_solve(F, projected_data);
//  dlog(1) << "projected_data=" << projected_data;

  // Faz o mapeamento para o vetor solucao 
  const DofMap & dof_map = sys.get_dof_map();
  std::vector<dof_id_type> dof_indices;
  dof_map.dof_indices (elem, dof_indices, var);

  // Cada processador seta os seus nos apenas
  dof_id_type f = sys.solution->first_local_index();
  dof_id_type l = sys.solution->last_local_index();
  for(uint i=0; i<n_dofs; i++)
  {
    dof_id_type dof_i = dof_indices[i];
    if ((f <= dof_i) && (dof_i < l )) 
      sys.solution->set(dof_i, projected_data(i));
  }
}

/**
 *
 *
 */
void ElemProjection::close_system() {
  SCOPELOG1(9);
  chimas_sync_check();
  sys.solution->close();
  sys.update();
}
