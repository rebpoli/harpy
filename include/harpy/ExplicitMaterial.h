#pragma once

#include "base/Global.h" 
#include "harpy/Material.h"

namespace libMesh { class MeshBase; class ExplicitSystem; class Elem; }

/**
 *   A class of material where the properties are projected explicitly 
 *   into the system.
 *
 *   This should not be instantiated directly, but inherited by 
 *   specific projecting classes that will define the variables etc.
 */
class ExplicitMaterial : public Material
{
public:
  ExplicitMaterial( suint sid_, const MaterialConfig & config_,
          ExplicitSystem & sys_ ) : Material( sid_, config_ ), system(sys_) {}

  ExplicitMaterial( Material & ref_material, ExplicitSystem & sys_ ) :
    Material( ref_material ), system(sys_) {}

  virtual void reinit( const Elem & elem_, uint side=255 ) ;
  virtual string hello() { return "ExplicitMaterial."; }

  void project( vector<double> & vals_qp, string vname="" );
  void project_tensor( vector<RealTensor> & vals_qp, string vname="" );

  void eval( vector<double> & vals_qp , string vname );

  void close_system();

protected:
  ExplicitSystem & system;
};
