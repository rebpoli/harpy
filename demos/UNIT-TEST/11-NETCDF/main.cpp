#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <exception>
#include <boost/mpi.hpp>
#include <netcdf.h>
#include <netcdf_par.h>
#include <libmesh/point.h>
#include <libmesh/tensor_value.h>

#include "harpy/Global.h"

using libMesh::Point;
using libMesh::RealTensorValue;

namespace mpi = boost::mpi;

// Variable types enum
enum class NC_VAR { TEMPERATURE, PRESSURE, VELOCITY, STRESS, DENSITY, VISCOSITY };

#define CHECK_NC(call) do { \
  int retval = call; \
  if (retval != NC_NOERR) flog << nc_strerror(retval); \
} while(0)

/**
 *
 */
class NetCDFWriter 
{
  private:
    mpi::communicator world;
    int ncid;

    // Dimension and variable IDs
    int time_dimid, point_dimid, vec3_dimid, ten9_dimid;
    int coord_varid, time_varid;

    // Variable maps
    map<NC_VAR, int> varids;

    // Static variable name mappings
    static map<string, NC_VAR> name_to_var_map;
    static map<NC_VAR, string> var_to_name_map;

  public:    

    uint n_times;
    uint n_points;

    // File state
    bool file_created;
    bool define_mode;

    explicit NetCDFWriter(uint n_times, uint n_points);
    ~NetCDFWriter();

    // File management
    void create_file(const string& filename);
    void close_file();

    // Variable definition methods
    void define_scalar_variable(const string & var_name, const string& units = "");
    void define_vector_variable(const string & var_name, const string& units = "");
    void define_tensor_variable(const string & var_name, const string& units = "");
    void finish_definitions();

    // Coordinate and time setting methods
    void set_coordinates_collective(const vector<Point>& points);
    void set_time_collective(const vector<double>& times);

    // Data setting methods
    void set_scalar(NC_VAR var, uint point_idx, uint time_idx, double value);
    void set_vector(NC_VAR var, uint point_idx, uint time_idx, const Point& vector);
    void set_tensor(NC_VAR var, uint point_idx, uint time_idx, const RealTensorValue& tensor);

    // Utility methods
    static NC_VAR get_var_from_name(const string& name);
    static string get_name_from_var(NC_VAR var);                 // debugging

  private:
    void setup_coordinate_variables();
    void add_global_attributes();
    void set_collective_access(int varid);
    void set_independent_access(int varid);
};

/*
 *
 *
 */

#include <ctime>
#include <iomanip>
#include <sstream>

// Initialize static maps
map<string, NC_VAR> NetCDFWriter::name_to_var_map = {
  {"Temperature", NC_VAR::TEMPERATURE},
  {"Pressure", NC_VAR::PRESSURE},
  {"Velocity", NC_VAR::VELOCITY},
  {"Stress", NC_VAR::STRESS},
  {"Density", NC_VAR::DENSITY},
  {"Viscosity", NC_VAR::VISCOSITY}
};

map<NC_VAR, string> NetCDFWriter::var_to_name_map = {
  {NC_VAR::TEMPERATURE, "Temperature"},
  {NC_VAR::PRESSURE, "Pressure"},
  {NC_VAR::VELOCITY, "Velocity"},
  {NC_VAR::STRESS, "Stress"},
  {NC_VAR::DENSITY, "Density"},
  {NC_VAR::VISCOSITY, "Viscosity"}
};

/** **/
NetCDFWriter::NetCDFWriter(uint n_times, uint n_points)
  : world(), n_times(n_times), n_points(n_points),
  file_created(false), define_mode(false) { }

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
  CHECK_NC(nc_def_dim(ncid, "ten9_comp", 9, &ten9_dimid));

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

/** **/
void NetCDFWriter::define_scalar_variable(const string & var_name, const string& units) 
{
  if (!define_mode) flog << "Not in define mode";

  int varid;
  int dimids[2] = {time_dimid, point_dimid};

  CHECK_NC(nc_def_var(ncid, var_name.c_str(), NC_DOUBLE, 2, dimids, &varid));

  // Add attributes
  if (!units.empty()) CHECK_NC(nc_put_att_text(ncid, varid, "units", units.length(), units.c_str()));

  auto var = get_var_from_name(var_name);
  varids[var] = varid;

  if (world.rank() == 0) ilog << "Defined scalar variable: " << var_name;
}

/** **/
void NetCDFWriter::define_vector_variable(const string & var_name, const string& units) 
{
  if (!define_mode) flog << "Not in define mode";

  int varid;
  int dimids[3] = {time_dimid, point_dimid, vec3_dimid};

  CHECK_NC(nc_def_var(ncid, var_name.c_str(), NC_DOUBLE, 3, dimids, &varid));

  // Add attributes
  if (!units.empty()) CHECK_NC(nc_put_att_text(ncid, varid, "units", units.length(), units.c_str()));

  auto var = get_var_from_name( var_name );
  varids[var] = varid;

  if (world.rank() == 0) ilog << "Defined vector variable: " << var_name << " (3 components)";
}

/** **/
void NetCDFWriter::define_tensor_variable( const string& var_name, const string& units ) 
{
  if (!define_mode) flog << "Not in define mode";

  int varid;
  int dimids[3] = {time_dimid, point_dimid, ten9_dimid};

  CHECK_NC(nc_def_var(ncid, var_name.c_str(), NC_DOUBLE, 3, dimids, &varid));

  // Add attributes
  if (!units.empty()) CHECK_NC(nc_put_att_text(ncid, varid, "units", units.length(), units.c_str()));

  auto var = get_var_from_name( var_name );
  varids[var] = varid;

  if (world.rank() == 0) ilog << "Defined tensor variable: " << var_name << " (3 components)";
}

/** **/
void NetCDFWriter::finish_definitions() 
{
  if (!define_mode) flog << "Not in define mode";

  CHECK_NC(nc_enddef(ncid));
  define_mode = false;

  if (world.rank() == 0) {
    ilog << "Finished variable definitions";
  }
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
void NetCDFWriter::set_scalar(NC_VAR var, uint point_idx, uint time_idx, double value) 
{
  if (define_mode) flog << "Still in define mode - call finish_definitions() first";

  auto it = varids.find(var);
  if (it == varids.end()) flog << "Scalar variable not defined: " << get_name_from_var(var);
  set_independent_access(it->second);

  luint start[2] = {time_idx, point_idx};
  luint count[2] = {1, 1};

  CHECK_NC(nc_put_vara_double(ncid, it->second, start, count, &value));
}

/** **/
void NetCDFWriter::set_vector(NC_VAR var, uint point_idx, uint time_idx, const Point& p) 
{
  if (define_mode) flog << "Still in define mode - call finish_definitions() first";

  luint start[3] = { time_idx, point_idx, 0};
  luint count[3] = {1, 1, 3};

  auto it = varids.find(var);
  if (it == varids.end()) flog << "Scalar variable not defined: " << get_name_from_var(var);
  auto vid = it->second;
  set_independent_access(vid);

  std::vector<double> vec = {p(0), p(1), p(2)};
  CHECK_NC( nc_put_vara_double( ncid, vid, start, count, vec.data() ) );
}

/** **/
void NetCDFWriter::set_tensor(NC_VAR var, uint point_idx, uint time_idx, const RealTensorValue& tensor) 
{
  if (define_mode) flog << "Still in define mode - call finish_definitions() first";

  luint start[3] = { time_idx, point_idx, 0};
  luint count[3] = {1, 1, 9};

  auto it = varids.find(var);
  if (it == varids.end()) flog << "Scalar variable not defined: " << get_name_from_var(var);
  auto vid = it->second;
  set_independent_access(vid);

  std::vector<double> flat(9); 
  for (uint i = 0; i < 3; ++i) 
  for (uint j = 0; j < 3; ++j) 
    flat[i*3 + j] = tensor(i,j);

  CHECK_NC( nc_put_vara_double( ncid, vid, start, count, flat.data() ) );
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
NC_VAR NetCDFWriter::get_var_from_name(const string& name) 
{
  auto it = name_to_var_map.find(name);
  if (it == name_to_var_map.end()) flog << "Unknown variable name: " + name;
  return it->second;
}
string NetCDFWriter::get_name_from_var(NC_VAR var) 
{
  auto it = var_to_name_map.find(var);
  if (it == var_to_name_map.end()) flog << "Unknown variable enum";
  return it->second;
}

#include <cmath>
#include <vector>

vector<double> generate_times(int n_times) {
  vector<double> times(n_times);
  for (int t = 0; t < n_times; ++t) {
    times[t] = t * 86400.0; // Daily intervals in seconds
  }
  return times;
}

double generate_temperature(int point_idx, int time_step) {
  double base_temp = 20.0;
  double spatial_var = 5.0 * sin(point_idx * 0.01);
  double temporal_var = 10.0 * sin(time_step * 0.1);
  return base_temp + spatial_var + temporal_var;
}

Point generate_velocity(int point_idx, int time_step) {
  double u = sin(point_idx * 0.02 + time_step * 0.05);
  double v = cos(point_idx * 0.03 + time_step * 0.07);
  double w = 0.1 * sin(point_idx * 0.01 + time_step * 0.03);
  return Point(u, v, w);
}

RealTensorValue generate_stress(int point_idx, int time_step) {
  RealTensorValue stress;
  double base = point_idx * 0.001 + time_step * 0.01;

  // Symmetric stress tensor
  stress(0,0) = 1000 + 100 * sin(base);        // XX
  stress(1,1) = 1000 + 80 * cos(base);         // YY
  stress(2,2) = 1000 + 60 * sin(base * 1.5);   // ZZ
  stress(0,1) = stress(1,0) = 50 * sin(base * 2);   // XY
  stress(0,2) = stress(2,0) = 30 * cos(base * 1.8); // XZ
  stress(1,2) = stress(2,1) = 20 * sin(base * 2.2); // YZ

  return stress;
}

int main(int argc, char** argv) {
    mpi::environment env(argc, argv);
    mpi::communicator world;

    // Parameters
    const int N_TIMES = 100;   // Reduced for example
    const int N_POINTS = 1000; // Reduced for example
    // Generate points
    vector<Point> points;
    for (uint i=0; i<N_POINTS; i++) points.emplace_back( i, i+1, i+2 );

    // Create NetCDF writer
    NetCDFWriter writer(N_TIMES, N_POINTS);

    if (world.rank() == 0) {
      ilog << "Starting NetCDF parallel write example";
      ilog << "Processes: " << world.size();
    }

    // Create file and define variables
    writer.create_file("example_data.nc");

    writer.define_scalar_variable("Temperature", "Kelvin");
    writer.define_scalar_variable("Pressure", "Pa");
    writer.define_vector_variable("Velocity", "m/s");
    writer.define_tensor_variable("Stress", "Pa");

    writer.finish_definitions();


    writer.set_coordinates_collective( points );

    // Generate and set time data
    auto times = generate_times(N_TIMES);
    writer.set_time_collective(times);

    // Main write loops (as requested in main function)
    if (world.rank() == 0) ilog << "Writing data...";

    // Each processor has its own points
    uint pts_per_rank = N_POINTS / world.size();  
    uint start_idx = world.rank() * pts_per_rank;
    uint end_idx   = start_idx + pts_per_rank; 

    std::cout << "Rank " << world.rank() << " handles points " 
      << start_idx << " .. " << end_idx-1 << "\n";



    // Write scalar data (temperature, pressure) independently
    for (uint time = 0; time < N_TIMES; ++time) 
    {
      if (world.rank() == 0 && (time % 10 == 0 || time < 5)) { ilog << "  Processing time step " << time + 1 << " of " << N_TIMES; }

      // Each processor writes data for its own local points
      for (uint i = start_idx; i < end_idx; ++i) 
      {
        double temp = generate_temperature(i, time);
        writer.set_scalar(NC_VAR::TEMPERATURE, i, time, temp);

        double pressure = 101325.0 + temp * 100; // Pressure depends on temperature
        writer.set_scalar(NC_VAR::PRESSURE, i, time, pressure);

        // Set vector values
        Point velocity = generate_velocity(i, time);
        writer.set_vector(NC_VAR::VELOCITY, i, time, velocity);

        // Set tensor values
        RealTensorValue stress = generate_stress(i, time);
        writer.set_tensor(NC_VAR::STRESS, i, time, stress);
      }
    }

    if (world.rank() == 0) {
      ilog << "\nData writing completed successfully!";
      ilog << "File structure:";
      ilog << "  Coordinates: X, Y, Z (collective)";
      ilog << "  Time: time_in_seconds (collective)";
      ilog << "  Scalars: temperature, pressure (independent)";
      ilog << "  Vectors: velocity_X, velocity_Y, velocity_Z (independent)";
      ilog << "  Tensors: stress_XX, stress_XY, ... stress_ZZ (independent)";
    }

//    // File automatically closed by destructor


    return 0;
}
