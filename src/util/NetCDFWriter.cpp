
#include "util/NetCDFWriter.h"

namespace util {

map<NC_VAR, DataInfo> PARAMS = {
    { TEMPERATURE, {SCALAR, "Temperature", "K",   "", -1} },
    { PRESSURE,    {SCALAR, "Pressure",    "Pa",  "", -1} },
    { VELOCITY,    {VEC3,   "Velocity",    "m/s", "", -1} },
    { STRESS,      {TEN9,    "Stress",      "Pa",  "", -1} }
};

/** **/
NetCDFWriter::NetCDFWriter(uint n_times, uint n_points)
  : world(), n_times(n_times), n_points(n_points),
  file_created(false), define_mode(false),
  vec3_dimid(0), ten9_dimid(0) { }

  /** **/
  NetCDFWriter::~NetCDFWriter() {
    if (file_created) close_file();
  }

/** **/
void NetCDFWriter::create_file(const string& filename) {
  if (file_created) flog << "File already created";

  MPI_Info info = MPI_INFO_NULL;

  // Create parallel NetCDF file
  CHECK_NC(nc_create_par(filename.c_str(), NC_NETCDF4 | NC_MPIIO,
        world, info, &ncid));

  // Define dimensions
  CHECK_NC(nc_def_dim(ncid, "time", n_times, &time_dimid));
  CHECK_NC(nc_def_dim(ncid, "point_idx", n_points, &point_dimid));
  CHECK_NC(nc_def_dim(ncid, "vec3_comp", 3, &vec3_dimid));

  setup_coordinate_variables();
  add_global_attributes();

  file_created = true;
  define_mode = true;

  if (world.rank() == 0) {
    ilog << "Created NetCDF file: " << filename;
    ilog << "Dimensions: time(" << n_times << "), point_idx(" << n_points << ")";
  }
}

/** **/
void NetCDFWriter::setup_coordinate_variables() 
{
  // Define coordinate variables
  CHECK_NC(nc_def_var(ncid, "time_in_seconds", NC_DOUBLE, 1, &time_dimid, &time_varid));

  int dimids[2] = {point_dimid, vec3_dimid};
  CHECK_NC(nc_def_var(ncid, "Coord", NC_DOUBLE, 2, dimids, &coord_varid));

  // Add coordinate attributes
  CHECK_NC(nc_put_att_text(ncid, time_varid, "units", 1, "s"));
  CHECK_NC(nc_put_att_text(ncid, coord_varid, "units", 1, "m"));
}

/** **/
void NetCDFWriter::add_global_attributes() 
{
  //    string title = "Parallel NetCDF Data";
  //    CHECK_NC(nc_put_att_text(ncid, NC_GLOBAL, "title",
  //                            title.length(), title.c_str()));

  //    auto creation_time = get_current_time_string();
  //    CHECK_NC(nc_put_att_text(ncid, NC_GLOBAL, "creation_date",
  //                            creation_time.length(), creation_time.c_str()));
}

/**
 *
 */
void NetCDFWriter::define_variable( NC_VAR var )
{
  if (!define_mode) flog << "Not in define mode";
  auto & di = get_di(var);

  vector<int> dimids;
  if ( di.type == SCALAR ) dimids = {time_dimid, point_dimid};
  if ( di.type == VEC3 )   dimids = {time_dimid, point_dimid, vec3_dimid};

  if ( di.type == TEN9 )   
  {
    if ( ! ten9_dimid )   // Create the dimension, if not existing
      CHECK_NC(nc_def_dim(ncid, "ten9_comp", 9, &ten9_dimid));

    dimids = {time_dimid, point_dimid, ten9_dimid};
  }

  ilog << "NAME:" << di.name;
  CHECK_NC(nc_def_var(ncid, di.name.c_str(), NC_DOUBLE, dimids.size(), dimids.data(), &di.vid));
  CHECK_NC(nc_put_att_text(ncid, di.vid, "units", di.unit.length(), di.unit.c_str()));


  // Useful annotations
  if ( di.type == TEN9 ) {
    string comp = "xx xy xz yx yy yz zx zy zz";
    CHECK_NC(nc_put_att_text(ncid, di.vid, "components", comp.size(), comp.c_str()));    
  }
  if ( di.type == VEC3 ) {
    string comp = "x y z";
    CHECK_NC(nc_put_att_text(ncid, di.vid, "components", comp.size(), comp.c_str()));    
  }

  if (world.rank() == 0) ilog << "Defined variable: " << di.name;
}

/** **/
void NetCDFWriter::finish_definitions() 
{
  if (!define_mode) flog << "Not in define mode";
  CHECK_NC(nc_enddef(ncid));
  define_mode = false;
  if (world.rank() == 0) ilog << "Finished variable definitions";
}

/** **/
void NetCDFWriter::set_coordinates_collective(const vector<Point>& points) 
{
  if (define_mode) flog << "Still in define mode - call finish_definitions() first";

  // Set collective access for coordinate variables
  set_collective_access(coord_varid);

  uint npoints = points.size();

  // Extract coordinates
  vector<double> flat( npoints*3 );

  for (uint i = 0; i < npoints; ++i) 
  for (uint j = 0; j < 3; ++j) 
    flat[i*3+j] = points[i](j);  

  // Write coordinates collectively
  luint start[2] = { 0 ,      0 };
  luint count[2] = { npoints, 3 };

  CHECK_NC(nc_put_vara_double(ncid, coord_varid, start, count, flat.data()));

  if (world.rank() == 0) ilog << "Set coordinates collectively";
}

/** **/
void NetCDFWriter::set_time_collective(const vector<double>& times) 
{
  if (define_mode) flog << "Still in define mode - call finish_definitions() first";

  if (times.size() != n_times) flog << "Times vector size doesn't match n_times";

  // Set collective access for time variable
  set_collective_access(time_varid);

  // Write time data collectively (all processes write same data)
  luint start[1] = {0};
  luint count[1] = {n_times};

  CHECK_NC(nc_put_vara_double(ncid, time_varid, start, count, times.data()));

  if (world.rank() == 0) ilog << "Set time collectively";
}

/** **/
void NetCDFWriter::set_value(NC_VAR var, double value) 
{
  if (define_mode) flog << "Still in define mode - call finish_definitions() first";
  auto & di = get_di(var);

  set_independent_access(di.vid);

  luint start[2] = { curr_time, curr_pt };
  luint count[2] = {1, 1};

  CHECK_NC(nc_put_vara_double(ncid, di.vid, start, count, &value));
}

/** **/
void NetCDFWriter::set_value(NC_VAR var, const Point& p) 
{
  if (define_mode) flog << "Still in define mode - call finish_definitions() first";
  auto & di = get_di(var);

  luint start[3] = { curr_time, curr_pt, 0};
  luint count[3] = {1, 1, 3};

  std::vector<double> vec = {p(0), p(1), p(2)};
  CHECK_NC( nc_put_vara_double( ncid, di.vid, start, count, vec.data() ) );
}

/** **/
void NetCDFWriter::set_value(NC_VAR var, const RealTensorValue& tensor) 
{
  if (define_mode) flog << "Still in define mode - call finish_definitions() first";
  auto & di = get_di(var);

  luint start[3] = { curr_time, curr_pt, 0};
  luint count[3] = {1, 1, 9};

  std::vector<double> flat(9); 
  for (uint i = 0; i < 3; ++i) 
  for (uint j = 0; j < 3; ++j) 
    flat[i*3 + j] = tensor(i,j);

  CHECK_NC( nc_put_vara_double( ncid, di.vid, start, count, flat.data() ) );
}

/** **/
void NetCDFWriter::set_collective_access(int varid) { CHECK_NC(nc_var_par_access(ncid, varid, NC_COLLECTIVE)); }
void NetCDFWriter::set_independent_access(int varid) { CHECK_NC(nc_var_par_access(ncid, varid, NC_INDEPENDENT)); }
void NetCDFWriter::close_file() 
{
  if (! file_created) return;
  CHECK_NC(nc_close(ncid));
  file_created = false;
  if (world.rank() == 0) ilog << "Closed NetCDF file";
}

} // ns
