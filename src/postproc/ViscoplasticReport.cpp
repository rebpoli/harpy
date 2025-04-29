
#include "postproc/ViscoplasticReport.h"
#include "solver/ViscoplasticSolver.h"
#include "libmesh/mesh.h"

#include "util/OutputOperators.h"

/**
 *
 */
void ViscoplasticReport::do_export()
{
  for ( Probe * probe : probes ) 
  {
    if ( probe->is_gauss() )
      export_by_face( dynamic_cast<GaussProbe &>(*probe) );
    else
      export_by_point( *probe );
  }
}

/**
 *
 */
void ViscoplasticReport::export_by_point( Probe & probe )
{
  SCOPELOG(1);

  dlog(1) << "Processing probe '" << probe.name << "' ...";
  MeshBase & mesh = solver.get_mesh();
  unique_ptr<PointLocatorBase> plocator = mesh.sub_point_locator();

  for ( Point & pt : probe.points ) 
  {
    const Elem * elem = (*plocator)(pt);
    if ( ! elem ) { wlog << "Element not found in " << Print(pt) << "."; continue; }

    dlog(1) << "Element found at " << Print(pt) << ": " << elem->id();

    ViscoPlasticMaterial * mat = solver.get_material( *elem );
    mat->value_at( pt, elem );

    
  }
}

/**
 *
 */
void ViscoplasticReport::export_by_face( GaussProbe & probe )
{
  SCOPELOG(1);
  dlog(1) << "Is gauss? " << probe.is_gauss();
  dlog(1) << probe;
}
