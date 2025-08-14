#pragma once

#include "base/Global.h"
#include "config/MaterialConfig.h"
#include <util/Serialize.h>

#include "libmesh/tensor_value.h"
#include "libmesh/vector_value.h"

using namespace libMesh;
struct ViscoplasticIFC;

struct VPProps 
{
  // Initialization
  void init_from_config( const MaterialConfig & config, const Point & pt );
  // Calculated properties
  inline double C_ijkl( uint i, uint j, uint k, uint l);
  // Calculate stresses from displacements
  void update( const RealVectorValue & U_, const RealTensor & GRAD_U_, double dt );
  // Updates the initial strain (typically after initialization solve)
  void update_initial_strain();
  
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
  RealTensor initial_strain;

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
  vector<double> von_mises, epskk, F, pressure;
};

/// Probe interface
struct ProbeIFC { 
  ProbeIFC() = default;
  ProbeIFC( const ProbeIFC & ) = delete;
  uint elem_id; Point pt; VPProps props; 
};

// Add serialization capabilities and MPI support
using ProbeByElemBaseMap = map< uint, vector<ProbeIFC *> >;
struct ProbeByElemMap : public ProbeByElemBaseMap
{
  ProbeByElemMap() : ProbeByElemBaseMap(), world(), rank(world.rank()), size(world.size()) {}

  // Main function that collects data from all processes and localizes to rank 0
  void localize_to_one(ProbeByElemMap & global_map) 
  {
    if (size == 1) { if (rank == 0) global_map = *this; return; }  // single proc.

    global_map.clear();
    if (rank != 0) _send(); else _receive( global_map );
    world.barrier();   // sync
  }
private:
  mpi::communicator world;
  int rank, size;
  friend class boost::serialization::access;

  // On rank != 0
  void _send() {  world.send(0, 0, *this);  return; }

  // On rank 0
  void _receive( ProbeByElemMap & global_map ) 
  {
    global_map = *this;  // add my own data

    for (int i = 1; i < size; ++i) 
    {
      ProbeByElemMap remote_map;
      world.recv(i, 0, remote_map);

      for (const auto& [ key, vec ] : remote_map) 
      {
        auto & tvec = global_map[key];  // create if non existing
        tvec.insert(tvec.end(), vec.begin(), vec.end());
      }
    }
  }

  // For boost
  template<class Archive>
  void serialize(Archive& ar, const unsigned int version) {
    UNUSED(version);
    ar & boost::serialization::base_object<ProbeByElemBaseMap>(*this);
  }
};

using ProbeByPnameByElemMap = map< string, ProbeByElemMap > ;

/** 
 *
 * Holds the data to feed viscoplastic properties into/from the material 
 *
 **/
struct ViscoplasticIFC
{


  ~ViscoplasticIFC();

  /// Data storage. By element, By Qp
  map< uint, vector<VPProps> > by_elem; 

  /// Dynamic, set after initialization
  vector<VPProps> * by_qp;

  /// 1 if the interface has been initialized
  bool valid;   

  void reinit( uint eid, uint nqp=0 );

  /// Getters
  VPProps & get( uint qp ) { return (*by_qp)[qp]; }

  /// Probes organzied by elemetns
  ProbeByPnameByElemMap probes_by_pname_by_elem;  // This map is MPI aware

  void add_probe_point( const MaterialConfig & config, string & name, uint eid, const Point & pt );

};

// Add serialization for the specific types
namespace boost { namespace serialization {
/** **/
template<class Archive>
void serialize(Archive & ar, VPProps & p, const unsigned int version)
{
  UNUSED(version);
  ar & p.lame_mu; ar & p.lame_lambda; ar & p.GRAD_U;
  ar & p.U; ar & p.sigtot; ar & p.sigeff;
  ar & p.temperature; ar & p.initial_temperature;
  ar & p.pressure; ar & p.initial_pressure;
  ar & p.von_mises; 
} 
/** **/
template<class Archive>
void serialize(Archive & ar, ProbeIFC & p, const unsigned int version)
{ UNUSED(version); ar & p.elem_id; ar & p.pt; ar & p.props;  } 

} }  // Namespaces


/** OUTPUT STREAMS **/
ostream& operator<<(ostream& os, const ViscoplasticIFC & m);
ostream& operator<<(ostream& os, const VPProps & m);

/** Inlines **/
inline double VPProps::C_ijkl( uint i, uint j, uint k, uint l) {
  const auto kd = Math::kronecker_delta;  // From Global

  return lame_lambda * ( kd(i,j)*kd(k,l) ) 
       + lame_mu * ( kd(i, k)*kd(j, l) + kd(i, l)*kd(j, k) );
}
