#pragma once

#include "base/Global.h"
#include "harpy/Material.h"

/**
 *
 */

namespace libMesh { class System; }

class ViscoPlasticMaterial : public Material
{
public:
  ViscoPlasticMaterial( suint sid_, const MaterialConfig & config, System & sys_ );
  virtual void init_fem();


  virtual void reinit( const Elem & elem );
  virtual void jacobian (const NumericVector<Number> & soln, SparseMatrix<Number> & jacobian );
  virtual void residual (const NumericVector<Number> & soln, NumericVector<Number> & residual );

private:
  void setup_variables();

  vector<vector<dof_id_type>> dof_indices_var;
  vector<vector< DenseSubMatrix<Number> >> Ke_var;
  vector< DenseSubVector<Number> > Fe_var;

};
