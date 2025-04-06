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

private:
  void setup_variables();

  vector<vector<dof_id_type>> dof_indices_var;
  vector<vector< DenseSubMatrix<Number> >> Ke_var;
  vector< DenseSubVector<Number> > Fe_var;

};
