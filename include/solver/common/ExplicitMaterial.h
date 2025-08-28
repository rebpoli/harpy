#pragma once

#include "harpy/Global.h" 
#include "config/MaterialConfig.h"

#include "libmesh/quadrature_gauss.h"
#include "libmesh/tensor_value.h"
#include "libmesh/fe_base.h"

namespace libMesh { class MeshBase; class ExplicitSystem; class Elem;  }

namespace solver { namespace viscoplastic { class ViscoPlasticMaterial ; } }


namespace solver {
namespace common {

using solver::viscoplastic::ViscoPlasticMaterial;
using config::MaterialConfig;

using namespace libMesh;

/**
 *   A class of material where the properties are projected explicitly 
 *   into the system.
 *
 *   This should not be instantiated directly, but inherited by 
 *   specific projecting classes that will define the variables etc.
 */
class ExplicitMaterial 
{
public:
//  ExplicitMaterial( suint sid_, const MaterialConfig * config_,
//          ExplicitSystem & sys_ ) : sid(sid_), config( config_ ), system(sys_) {}

  ExplicitMaterial( ViscoPlasticMaterial * ref_material, ExplicitSystem & sys_ );
  virtual ~ExplicitMaterial() {}

  virtual void reinit( const Elem & elem_, uint side=255 ) ;
  virtual string hello() { return "ExplicitMaterial."; }

  void project( vector<double> & vals_qp, string vname="" );
  void project_tensor( vector<RealTensor> & vals_qp, string vname="" );
  void project_tensor_invariants( vector<RealTensor> & vals_qp, string vname="" );

  void eval( vector<double> & vals_qp , string vname );
  double eval( const Point & pt, string vname );

  void close_system();

  // Members
  ViscoPlasticMaterial * refmat;
  const MaterialConfig * config;
  ExplicitSystem & system;

  QGauss qrule;
  const Elem * elem;      /// The element that the material has been reinit'ed to
  unique_ptr<FEBase> fe;  /// The finite element object to hold shape funtions, jxw, etc

  string name;
  suint sid;   
};

/**
 *
 */
class ExplicitThermalMaterial : public ExplicitMaterial
{
public:
  ExplicitThermalMaterial( ViscoPlasticMaterial * ref_material, ExplicitSystem & sys_ ) :
    ExplicitMaterial( ref_material, sys_ ) {}

  virtual void init_fem();
};

/**
 *
 */
class ExplicitPressureMaterial : public ExplicitMaterial
{
public:
  ExplicitPressureMaterial( ViscoPlasticMaterial * ref_material, ExplicitSystem & sys_ ) :
    ExplicitMaterial( ref_material, sys_ ) {}
  virtual void init_fem();
};

}} // ns
