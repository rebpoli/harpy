
#include "postproc/ViscoplasticReport.h"
#include "solver/ViscoplasticSolver.h"
#include "material/ViscoPlasticMaterial.h"
#include "harpy/Timestep.h"
#include "util/OutputOperators.h"

#include "libmesh/mesh.h"


/**
 *
 */
ViscoplasticReport::ViscoplasticReport( ViscoplasticSolver & solver_ ) : solver(solver_) 
{ }

/**
 *    Register probes intothe material interface.
 */
void ViscoplasticReport::init()
{
  for ( Probe * probe : probes ) 
  {
    CsvFile1 ofile(probe->filename, "\t", false);
    ofile << "Time" << "Timestep" << "Var" << "X" << "Y" << "Z" << "Value" << endrow;

    init_material( *probe );
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
}

/**
 *
 */
void ViscoplasticReport::export_by_point( Probe & probe )
{
//  SCOPELOG(1);

//  dlog(1) << "Processing probe '" << probe.name << "' ...";
//  CsvFile1 ofile(probe.filename);

//  double time = solver.ts.time;
//  double t_step = solver.ts.t_step;

//  // Export the probe points stored in the material interfaces
//  for ( auto & [ sid, _ ] : solver.material_by_sid )
//  {
//    ViscoPlasticMaterial * mat = solver.get_material(sid);
//    if ( ! mat->vp_ifc.probes_by_pname.count(probe.name) ) continue;

//    auto & vec = mat->vp_ifc.probes_by_pname.at( probe.name ) ;
//    for ( auto & probe_ifc : vec ) 
//    {
//      auto & pt = probe_ifc.pt;
//      auto & p = probe_ifc.props;

//      auto & U = p.U;
//      ofile << t_step << time << "UX" << pt(0) << pt(1) << pt(2) << U(0) << endrow;
//      ofile << t_step << time << "UY" << pt(0) << pt(1) << pt(2) << U(1) << endrow;
//      ofile << t_step << time << "UZ" << pt(0) << pt(1) << pt(2) << U(2) << endrow;

//      auto & stot = p.sigtot;
//      ofile << t_step << time << "STOTXX" << pt(0) << pt(1) << pt(2) << stot(0,0) << endrow;
//      ofile << t_step << time << "STOTYY" << pt(0) << pt(1) << pt(2) << stot(1,1) << endrow;
//      ofile << t_step << time << "STOTZZ" << pt(0) << pt(1) << pt(2) << stot(2,2) << endrow;
//      ofile << t_step << time << "STOTYZ" << pt(0) << pt(1) << pt(2) << stot(1,2) << endrow;
//      ofile << t_step << time << "STOTXZ" << pt(0) << pt(1) << pt(2) << stot(0,2) << endrow;
//      ofile << t_step << time << "STOTXY" << pt(0) << pt(1) << pt(2) << stot(0,1) << endrow;

//      auto & T = p.temperature;
//      ofile << t_step << time << "T" << pt(0) << pt(1) << pt(2) << T << endrow;
//    }
//  }
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
