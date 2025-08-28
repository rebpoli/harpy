#pragma once

#include "harpy/Global.h"
#include "config/MaterialConfig.h"
#include <util/Serialize.h>
#include "libmesh/tensor_value.h"
#include "libmesh/vector_value.h"

namespace solver {
namespace viscoplastic {

using config::MaterialConfig;
using config::CreepMD;
using libMesh::Point;
using libMesh::RealVectorValue;
using libMesh::RealTensor;

/**
 *   Full set of data of a Viscoplastic Intragation Point
 */
struct VPProps 
{
  // Initialization
  void init_from_config( const MaterialConfig & config, const Point & pt );
  // Calculated properties
  inline double C_ijkl( uint i, uint j, uint k, uint l);
  // Calculate stresses from displacements
  void update( const RealVectorValue & U_, const RealTensor & GRAD_U_, double dt );

  // Static variables
  double lame_mu, lame_lambda, alpha_d, beta_e;

  // Poroelasticity
  double biot, density;

  CreepMD creep_md; // Defined in MaterialConfig

  // State variables
  RealTensor GRAD_U;
  RealVectorValue U;
  double temperature, initial_temperature;
  double pressure, initial_pressure;

  // Stresses
  RealTensor sigtot, sigeff, deviatoric;
  double von_mises, epskk, F;

  // Stress Initialization
  RealTensor initial_stress;

  // Plasticity
  RealTensor plastic_strain, plastic_strain_n, plastic_strain_k;
  RealTensor plastic_strain_rate;

  // Debug
  double grad_norm;
  double sigeff0, plast0;
  RealTensor plast0_t;
};


/// Transposed version of the properties
struct PropsTranspose
{
  PropsTranspose( vector<VPProps> * by_qp ) ;

  vector<RealTensor> sigtot, sigeff, sigeff_terz;
  vector<RealTensor> deviatoric, plastic_strain, plastic_strain_rate;
  vector<RealTensor> initial_stress;
  vector<double> von_mises, epskk, F, pressure;
};

/** 
 * Inlines 
 **/
inline double VPProps::C_ijkl( uint i, uint j, uint k, uint l) {
  const auto kd = Math::kronecker_delta;  // From Global

  return lame_lambda * ( kd(i,j)*kd(k,l) ) 
       + lame_mu * ( kd(i, k)*kd(j, l) + kd(i, l)*kd(j, k) );
}

/** 
 * Output streams 
 **/
ostream& operator<<(ostream& os, const VPProps & m);

}} // ns

/**
 * Add serialization capability for VPProps
 */
namespace boost {
namespace serialization {
  using solver::viscoplastic::VPProps;

  /** **/
  template<class Archive>
  void serialize(Archive & ar, VPProps & p, const unsigned int /*ver*/)
  {
    ar & p.lame_mu; ar & p.lame_lambda; ar & p.GRAD_U;
    ar & p.U; ar & p.sigtot; ar & p.sigeff;
    ar & p.temperature; ar & p.initial_temperature;
    ar & p.pressure; ar & p.initial_pressure;
    ar & p.von_mises; 
    ar & p.initial_stress; 
    ar & p.plastic_strain; 
  } 
} }  // Namespaces

