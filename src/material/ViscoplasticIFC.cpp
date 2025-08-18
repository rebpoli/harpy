
#include "material/ViscoplasticIFC.h"
#include "util/OutputOperators.h"

ViscoplasticIFC::~ViscoplasticIFC()
{
  for ( auto & [ _, m1 ] : probes_by_pname_by_elem )
  for ( auto & [ _, vec ] : m1 )
  for ( auto p : vec ) delete(p);
}

/**
 *
 */
void ViscoplasticIFC::reinit( uint eid, uint nqp )
{
  by_qp = &( by_elem[eid] );

  if ( nqp )
  if ( ! by_qp->size() ) 
  {
    by_qp->resize(nqp); 
    valid=0; 
  }
}

/**
 *
 */
void ViscoplasticIFC::add_probe_point( const MaterialConfig & config, string & name, uint eid, const Point & pt )
{
  VPProps props;
  props.init_from_config (config, pt );

  ProbeIFC * probe_ifc = new ProbeIFC() ;
  probe_ifc->elem_id = eid;
  probe_ifc->pt = pt;
  probe_ifc->props = props;

  auto & m1 = probes_by_pname_by_elem[name];
  auto & vec_e = m1[eid];
  vec_e.push_back( probe_ifc );   ///PROBE
};

/**
 *
 */
InitVPPropsByElemMap ViscoplasticIFC::snapshot_initial_strain()
{
  InitVPPropsByElemMap ret;
  for ( auto & [k,vec] : by_elem )
  {
    auto & dst = ret[k];
    for ( auto & v : vec ) 
    {
      // Create the snapshot on the fly
      InitVPProps e;
      e.initial_strain = v.initial_strain;
      // Register the snapshot in the map
      dst.push_back( e );
    }
  }
  return ret;
};

/**
 *
 */
void ViscoplasticIFC::save_initial_strain( const string & filename )
{
  // This should be done only in rank=0
  InitVPPropsByElemMap ivmap = snapshot_initial_strain();
  ivmap.save(filename);
}

void ViscoplasticIFC::load_initial_strain( const string & filename )
{
  InitVPPropsByElemMap ivmap;
  ivmap.load(filename);

  // TODO : add some validations here
  for ( auto & [k,vec] : by_elem )
  {
    auto & src = ivmap[k];
    for ( uint i=0; i<vec.size(); i++ )
      vec[i].initial_strain = src[i].initial_strain;
  }
}

/**
 *
 */
ostream& operator<<(ostream& os, const ViscoplasticIFC & m)
{
  os << "VISCOPLASTIC INTERFACE:" << endl;

  for ( auto & [ eid, pvec ] : m.by_elem )
  {
    os << "    EID:" << eid;
    os << "          [" << endl;
    for ( auto & p : pvec ) {
      os << "             lame_mu:" << setw(15) << p.lame_mu << endl;
      os << "         lame_lambda:" << setw(15) << p.lame_lambda << endl;
    }
    os << "          ]" << endl;
  }

  return os;
}

