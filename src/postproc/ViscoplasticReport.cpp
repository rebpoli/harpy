
#include "postproc/ViscoplasticReport.h"
#include "solver/ViscoplasticSolver.h"
#include "material/ViscoPlasticMaterial.h"
#include "harpy/Timestep.h"
#include "util/OutputOperators.h"
#include "util/String.h"

#include "libmesh/mesh.h"


/**
 *
 */
ViscoplasticReport::ViscoplasticReport( ViscoplasticSolver & solver_ ) :
  solver(solver_) , scalars_fn("run/csv/scalars.csv")
{ }

/**
 *    Register probes intothe material interface.
 */
void ViscoplasticReport::init()
{
  /// PROBE CSV FILES
  for ( Probe * probe : probes ) 
  {
    CsvFile1 ofile(probe->filename, "\t", false);
    ofile << "Time" << "Timestep" << "Var" << "X" << "Y" << "Z" << "Value" << endrow;
    init_material( *probe );
  }

  /// SCALAR CSV FILE
  {
    CsvFile1 ofile(scalars_fn, "\t", false);
    ofile << "Time" << "Timestep" << "Var" << "Value" << endrow;
  }
}

/**
 *
 */
void ViscoplasticReport::init_material( Probe & probe  )
{
  SCOPELOG(1);
  if ( probe.is_gauss() ) return;  

  MeshBase & mesh = solver.get_mesh();
  unique_ptr<PointLocatorBase> plocator = mesh.sub_point_locator();
  plocator->enable_out_of_mesh_mode();

  for ( uint pt_i=0; pt_i<probe.points.size(); pt_i++ ) 
  {
    Point & pt = probe.points[pt_i];
    const Elem * elem = (*plocator)(pt);
    if ( ! elem ) { wlog << "Element not found in " << Print(pt) << "."; continue; }

    // The probes only exist in the element processor
    if ( ! elem->active() ) continue;
    if ( elem->processor_id() != mesh.processor_id() ) continue;

    dlog(1) << "> Element found at " << Print(pt) << ": " << elem->id();

    ViscoPlasticMaterial * mat = solver.get_material( *elem );
    ViscoplasticIFC & vp_ifc = mat->vp_ifc;
    vp_ifc.add_probe_point( mat->config, probe.name, elem->id(), pt );
  }
}

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

  export_scalars();
}

/**
 *   Export scalars from the ViscoPlastic System
 */
void ViscoplasticReport::export_scalars()
{
  double time = solver.ts.time;
  double t_step = solver.ts.t_step;
  auto & system = solver.system;

  CsvFile1 ofile(scalars_fn);

  for ( auto & sv : MODEL->boundary_config.scalars )
  {
    uint vid = system.variable_number(sv.name);
    vector<dof_id_type> dofi;
    system.get_dof_map().SCALAR_dof_indices (dofi, vid);
    if ( dofi.size() > 1 ) flog << "Scalar variable should have only one DOF.";

    // I assume that the scalars are in all processors ...
    double val = (*system.current_local_solution)(dofi[0]);
    using namespace harpy_string;

    ofile << t_step << time << sv.name << CSVSci(val) << endrow;
  }
}
  
/**
 *
 */
void ViscoplasticReport::export_by_point( Probe & probe )
{
  SCOPELOG(1);

  dlog(1) << "Processing probe '" << probe.name << "' ...";
  CsvFile1 ofile(probe.filename);

  double time = solver.ts.time;
  double t_step = solver.ts.t_step;

  // Export the probe points stored in the material interfaces
  for ( auto & [ sid, _ ] : solver.material_by_sid )
  {
    ViscoPlasticMaterial * mat = solver.get_material(sid);

    ProbeByElemMap & local_map = mat->vp_ifc.probes_by_pname_by_elem[probe.name];

    // It only returns something in processor rank=0 (root)
    ProbeByElemMap global_map;
    local_map.localize_to_one( global_map );

    for ( auto & [ eid, vec ] : global_map ) 
    for ( auto & probe_ifc : vec ) 
    {
      auto & pt = probe_ifc->pt;
      auto & p = probe_ifc->props;

      auto & U = p.U;
      ofile << t_step << time << "UX" << pt(0) << pt(1) << pt(2) << U(0) << endrow;
      ofile << t_step << time << "UY" << pt(0) << pt(1) << pt(2) << U(1) << endrow;
      ofile << t_step << time << "UZ" << pt(0) << pt(1) << pt(2) << U(2) << endrow;

      auto & stot = p.sigtot;
      ofile << t_step << time << "STOTXX" << pt(0) << pt(1) << pt(2) << stot(0,0) << endrow;
      ofile << t_step << time << "STOTYY" << pt(0) << pt(1) << pt(2) << stot(1,1) << endrow;
      ofile << t_step << time << "STOTZZ" << pt(0) << pt(1) << pt(2) << stot(2,2) << endrow;
      ofile << t_step << time << "STOTYZ" << pt(0) << pt(1) << pt(2) << stot(1,2) << endrow;
      ofile << t_step << time << "STOTXZ" << pt(0) << pt(1) << pt(2) << stot(0,2) << endrow;
      ofile << t_step << time << "STOTXY" << pt(0) << pt(1) << pt(2) << stot(0,1) << endrow;

      auto & T = p.temperature;
      ofile << t_step << time << "T" << pt(0) << pt(1) << pt(2) << T << endrow;
      ofile << t_step << time << "INITIAL_T" << pt(0) << pt(1) << pt(2) << p.initial_temperature << endrow;
    }
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
