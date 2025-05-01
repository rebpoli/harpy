#pragma once

#include "base/Global.h"
#include "libmesh/tensor_value.h"
#include "libmesh/vector_value.h"
#include "config/MaterialConfig.h"

using namespace libMesh;

/** 
 *
 * Holds the data to feed viscoplastic properties into/from the material 
 *
 **/
struct ViscoplasticIFC
{

  struct Props 
  {
    // Initialization
    void init_from_config( const MaterialConfig & config, const Point & pt );
    // Calculated properties
    inline double C_ijkl( uint i, uint j, uint k, uint l);

    // Static variables
    double lame_mu, lame_lambda, alpha_d, beta_e;
    double creep_carter_a, creep_carter_q, creep_carter_n;

    // State variables
    RealTensor GRAD_U;
    RealVectorValue U;
    double temperature;

    // Stresses
    RealTensor sigtot, sigeff, deviatoric;
    double von_mises, epskk;

    // Plasticity
    RealTensor plastic_strain, plastic_strain_n;
    RealTensor plastic_strain_rate;
  };

  /// Transposed version of the properties
  struct PropsTranspose
  {
    PropsTranspose( vector<Props> * by_qp ) ;

    vector<RealTensor> sigtot, sigeff, deviatoric, plastic_strain, plastic_strain_rate;
    vector<double> von_mises, epskk;
  };

  ~ViscoplasticIFC();

  /// Data storage. By element, By Qp
  map< uint, vector<Props> > by_elem; 

  /// Dynamic, set after initialization
  vector<Props> * by_qp;

  /// 1 if the interface has been initialized
  bool valid;   

  void reinit( uint eid, uint nqp=0 );

  /// Getters
  Props & get( uint qp ) { return (*by_qp)[qp]; }

  /// Probe interface
  struct ProbeIFC { 
    ProbeIFC() = default;
    ProbeIFC( const ProbeIFC & ) = delete;
    uint elem_id; Point pt; Props props; 
  };
  map< string, vector<ProbeIFC *> > probes_by_pname;

  /// Probes organzied by elemetns
  map< uint, vector<ProbeIFC *> > probes_by_elem; 

  void add_probe_point( const MaterialConfig & config, string & name, uint eid, const Point & pt );

};

/** OUTPUT STREAMS **/
ostream& operator<<(ostream& os, const ViscoplasticIFC & m);
ostream& operator<<(ostream& os, const ViscoplasticIFC::Props & m);

/** Inlines **/
inline double ViscoplasticIFC::Props::C_ijkl( uint i, uint j, uint k, uint l) {
  const auto kd = Math::kronecker_delta;  // From Global

  return lame_lambda * ( kd(i,j)*kd(k,l) ) 
       + lame_mu * ( kd(i, k)*kd(j, l) + kd(i, l)*kd(j, k) );
}
