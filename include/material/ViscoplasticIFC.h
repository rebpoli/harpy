#pragma once

#include "base/Global.h"
#include "libmesh/tensor_value.h"

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
    double lame_mu, lame_lambda, alpha_d, beta_e;
    double creep_carter_a, creep_carter_q, creep_carter_n;

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

  /// Data storage. By element, By Qp
  map< uint, vector<Props> > by_elem; 
                                       
  /// Dynamic, set after initialization
  vector<Props> * by_qp;

  /// 1 if the interface has been initialized
  bool valid;   

  void reinit( uint eid, uint nqp );

  /// Getters
  Props & get( uint qp ) { return (*by_qp)[qp]; }
};

/** OUTPUT STREAMS **/
ostream& operator<<(ostream& os, const ViscoplasticIFC & m);
ostream& operator<<(ostream& os, const ViscoplasticIFC::Props & m);
