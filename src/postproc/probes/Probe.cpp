
#include "postproc/probes/Probe.h"
#include "util/String.h"
#include "util/OutputOperators.h"

#include "libmesh/string_to_enum.h"
#include "libmesh/enum_to_string.h"

using namespace util;

namespace postproc {
namespace probes {

using config::InoutConfig;
using util::operator<<;

/**
 *
 *
 */
ProbeCol::ProbeCol() 
{
  InoutConfig config;
  for ( auto & pconf : config.probes )
    Probe::Factory( *this, pconf  );
}

/**
 *
 */
Probe * Probe::Factory( vector<Probe *> & ret, ProbeConfig & config )
{
  if ( ! config.type ) flog << "Undefined type for probe '" << config.name << "'.";

  string & type = *(config.type);

  if ( iequals( type, "radial" ) ) 
    return ret.emplace_back( new RadialProbe(config) );

  if ( iequals( type, "planar" ) ) 
    return ret.emplace_back( new PlanarProbe(config) );

  if ( iequals( type, "linear" ) ) 
    return ret.emplace_back( new LinearProbe(config) );

  if ( iequals( type, "gauss"  ) ) 
    return ret.emplace_back( new GaussProbe(config) );

  flog << "Unknown probe typp '" << type << "'.";
  return 0;
}

/**
 *
 *
 */
Probe::Probe( ProbeConfig & config ) : 
  name(config.name),
  filename("run/csv/" + name + ".csv"),
  filename_pq("run/csv/" + name + "-pq.csv")
{ }

/**
 *   Create a planar probe from the config.
 */
PlanarProbe::PlanarProbe( ProbeConfig & config ) : Probe(config)
{
  if ( ! config.p0 )   flog << "P0 not defined for probe '" << config.name << "'!";
  if ( ! config.p1 )   flog << "P1 not defined for probe '" << config.name << "'!";
  if ( ! config.p2 )   flog << "P2 not defined for probe '" << config.name << "'!";
  if ( ! config.npts ) flog << "NPTS not defined for probe '" << config.name << "'!";

  Point p0 = *config.p0,
        p1 = *config.p1,
        p2 = *config.p2;
  
  // Grab nx and ny (minimum of 5 in each direction)
  uint base = max( 5u, uint( sqrt( double(*config.npts) ) ) ) ,
       nx = base ,
       ny = base ;  

  Point e1 = p1 - p0;           // edge 1
  Point raw_e2 = p2 - p0;       // edge 2 (to be made independent of e1)

  // Remove component of raw_e2 along e1 for a clean sweep basis
  // operator* is a dot product
  double alpha = (raw_e2 * e1) / (e1 * e1);
  Point e2 = raw_e2 - e1 * alpha;

  // Increments
  Point  step1  = e1 / double(nx) ;
  Point  step2  = e2 / double(ny) ;

  Point row_start = p0;
  for (uint j = 0; j <= ny; ++j) {
      Point p = row_start;
      for (uint i = 0; i <= nx; ++i) {
          points.emplace_back(p);
          p += step1;
      }
      row_start += step2;
  }

  dlog(1) << "PLANAR PROBE '" << config.name << "': " << *this;
}

/**
 *
 * Cria o probe radial com uma lista de pontos.
 *
 */
RadialProbe::RadialProbe( ProbeConfig & config ) : Probe(config)
{
  if ( ! config.dtheta ) flog << "dtheta not defined for probe '" << config.name << "'!";
  if ( ! config.radius ) flog << "radius not defined for probe '" << config.name << "'!";
  if ( ! config.center ) flog << "center not defined for probe '" << config.name << "'!";
  if ( ! config.normal ) flog << "normal not defined for probe '" << config.name << "'!";

  double dth = *config.dtheta / 180 * PI;
  double r = *config.radius;
  Point c = *config.center;
  Point n = *config.normal;
  n /= n.norm(); 

  Point a = (abs(n(2)) < 0.99) ? Point(0, 0, 1) : Point(0, 1, 0);

  Point u = a.cross(n);
  u /= u.norm();

  Point v = n.cross(u); 

  for ( double th=0; th<2*PI; th+=dth)
  {
    Point p = c + 
              r * std::cos(th) * u + 
              r * std::sin(th) * v;

    points.push_back( p );
  }
}

/**
 *
 */
LinearProbe::LinearProbe( ProbeConfig & config ) : Probe(config)
{
  if ( ! config.to )   flog << "to not defined for probe '" << config.name << "'!";
  if ( ! config.from ) flog << "from not defined for probe '" << config.name << "'!";
  if ( ! config.npts ) flog << "npts not defined for probe '" << config.name << "'!";

  Point to = *config.to;
  Point from = *config.from;
  uint n = *config.npts;

  Point dp = ( to - from ) / double(n-1);
  Point curr = from;
  for ( uint i = 0 ; i < n ; i++ ) {
    points.push_back( curr );
    curr = curr + dp;
  }
}

/**
 *
 */
GaussProbe::GaussProbe( ProbeConfig & config ):
  Probe(config), order(static_cast<Order>(FIRST))
{ 
  if ( ! config.order )       flog << "order not defined for probe '" << config.name << "'!";
  if ( ! config.boundaries.size() )  flog << "No boundary defined for probe '" << config.name << "'!";

  order = Utility::string_to_enum<Order>( *config.order ) ;
  boundaries = config.boundaries;
}

///**
// *
// *
// */
//void GaussProbe::eval( System & sys, vector<string> vars ) 
//{
//  dlog1(1) << "Evaluating GaussProbe to file '"<< filename <<"' ... ORDER: "<<order;
//  CsvFile1 ofile(filename);

//  const MeshBase & mesh = sys.get_mesh();
//  const BoundaryInfo & bi = mesh.get_boundary_info();
//  auto target_bid = bi.get_id_by_name( bname );
//  vector<double> soln;
//  sys.solution->localize(soln);
//  const DofMap & dof_map = sys.get_dof_map();

//  /// Para cada variavel ...
//  vector<uint> ivars = get_ivars(sys, vars);
//  for ( uint vari : ivars ) {
//    /// Inicializa FE
//    FEType fetype = sys.variable_type( vari );
//    unique_ptr<FEBase> fe = FEBase::build( 3, fetype );
//    QGauss qgauss(2, order);
//    fe->attach_quadrature_rule ( &qgauss );
//    const vector<vector<Real> > & phi = fe->get_phi();
//    const vector<Point> & xyz = fe->get_xyz(); 

//    /// Para cada elemento, procura a borda
//    for (const auto & elem : mesh.active_element_ptr_range()) 
////    if ( SD_CONT.count( elem->subdomain_id() ) ) // só elementos do continuo
//    for (auto side : elem->side_index_range()) {
//      vector<boundary_id_type> vec;
//      bi.boundary_ids(elem, side, vec );
//      for ( auto bid : vec ) {
//        if ( bid != target_bid ) continue;
//        vector<dof_id_type> dofi;
//        dof_map.dof_indices (elem, dofi, vari);

//        if ( ! dofi.size() ) continue; // elemento nao tem graus de liberdade nessa variavel

//        fe->reinit(elem, side);
//        for (uint qp=0; qp<qgauss.n_points(); qp++) {
//          double v=0;
//          for (uint l=0; l<dofi.size(); l++) 
//            v += phi[l][qp] * soln[ dofi[l] ];

//          const Point & pt = xyz[qp];

//          ofile << time << sys.variable_name(vari) << pt(0) << pt(1) << pt(2) << v << endrow;
//        }
//      }
//    }
//  }
//}

///**
// *
// *
// */
//void GaussProbe::eval_normal(System & sys, string src_var, string trg_var, double scale ) 
//{
//  dlog1(1) << "Evaluating GaussProbe to file '"<< filename <<"' ...";
//  CsvFile1 ofile(filename);

//  const MeshBase & mesh = sys.get_mesh();
//  const BoundaryInfo & bi = mesh.get_boundary_info();
//  auto target_bid = bi.get_id_by_name( bname );
//  vector<double> soln;
//  sys.solution->localize(soln);
//  const DofMap & dof_map = sys.get_dof_map();
//  uint vari = sys.variable_number( src_var );
//  /// Inicializa FE
//  FEType fetype = sys.variable_type( vari );
//  unique_ptr<FEBase> fe = FEBase::build( 3, fetype );
//  QGauss qgauss(2, order);
//  fe->attach_quadrature_rule ( &qgauss );
//  const vector<vector<Real> > & phi = fe->get_phi();
//  const vector<vector<RealGradient> > & dphi = fe->get_dphi();
//  const vector<Point> & xyz = fe->get_xyz(); 
//  const vector<Point> & normals = fe->get_normals(); 

//  /// Para cada elemento, procura a borda
//  for (const auto & elem : mesh.active_element_ptr_range()) 
//  for (auto side : elem->side_index_range()) {
//    vector<boundary_id_type> vec;
//    bi.boundary_ids(elem, side, vec );
//    for ( auto bid : vec ) {
//      if ( bid != target_bid ) continue;
//      vector<dof_id_type> dofi;
//      dof_map.dof_indices (elem, dofi, vari);
//      if ( ! dofi.size() ) continue; // elemento nao tem graus de liberdade nessa variavel

//      fe->reinit(elem, side);
//      for (uint qp=0; qp<qgauss.n_points(); qp++) {
//        double v=0;
//        Gradient dv = 0;
//        for (uint l=0; l<dofi.size(); l++) {
//          v += phi[l][qp] * soln[ dofi[l] ];
//          dv.add_scaled ( dphi[l][qp], soln[ dofi[l] ] );
//        }

//        const Point & pt = xyz[qp];
//        double dvdn = dv * normals[qp] * scale;
//        ofile << time << trg_var << pt(0) << pt(1) << pt(2) << dvdn << endrow;
//      }
//    }
//  }
//}

/**
 *
 *
 */
void PlanarProbe::print( ostream& os ) const
{ os << setw(20) << "PlanarProbe:" << setw(20) << "" << Print(points) << endl; }
//
void RadialProbe::print( ostream& os ) const
{ os << setw(20) << "RadialProbe:" << setw(20) << "" << Print(points) << endl; }
//
void LinearProbe::print( ostream& os ) const
{ os << setw(20) << "LinearProbe:" << setw(20) << "" << Print(points) << endl; }
//
void GaussProbe::print( ostream& os ) const
{ 
  os << setw(20) << "GaussProbe:" << endl;
  os << setw(30) << "Order" << setw(10) << "" <<  Utility::enum_to_string(order) << endl;
  os << setw(30) << "Boundaries" << setw(10) << "" << boundaries << endl;
}
//
ostream& operator<<(ostream& os, const ProbeCol & m)
{
  os << "Probe Collection (ProbeCol):"<< endl;
  for ( auto & p : m ) os << *p ; 
  return os;
}
//
ostream& operator<<(ostream& os, const Probe & m) 
{
  m.print(os);
  return os;
}

//ostream& operator<<(ostream& os, const vector<Probe *> & m) {
//  for ( const Probe * probe : m ) {
//    os << "Probe '"<<probe->get_name()<<"':" << endl;
//    os << "    " << *probe;
//  }
//  return os;
//}

}} // ns
