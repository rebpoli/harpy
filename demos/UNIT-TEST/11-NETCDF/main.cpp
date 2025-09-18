#include <netcdf.h>
#include <vector>
#include <string>
#include <unordered_map>
#include <stdexcept>
#include <iostream>
#include "harpy/HarpyInit.h"
#include "libmesh/point.h"
#include "libmesh/tensor_value.h"
#include "harpy/Global.h"
#include "util/Stopwatch.h"

using namespace std;
using namespace libMesh;
using namespace harpy;
using namespace util;

/**
 *
 */
class NetCDFWriter 
{
  public:
    NetCDFWriter(const string &filename);
    ~NetCDFWriter() ;
    void define_dimensions(size_t npoints) ;
    void define_scalar(const string &name) ;
    void define_vector(const string &name, int ncomp) ;
    void define_tensor(const string &name) ;
    void end_define() ;
    void set_time(size_t t) ;
    void set_coord(size_t point_idx, double x, double y, double z) ;
    void set_coord(size_t point_idx, const Point &p) ;
    void set_scalar(const string &name, size_t point_idx, double val) ;
    void set_vector(const string &name, size_t point_idx, const vector<double> &vals) ;
    void set_point(const string &name, size_t point_idx, const vector<double> &vals) ;
    void set_point(const string &name, size_t point_idx, const Point &p) ;
    void set_tensor(const string &name, size_t point_idx, const vector<vector<double>> &vals) ;
    void set_tensor(const string &name, size_t point_idx, const RealTensor &T) ;
    void define_unit(const string &var_name, const string &units) ;

  private:
    int ncid;
    int dim_time, dim_points;
    int var_x, var_y, var_z;
    size_t current_time;
    bool data_mode = false;

    unordered_map<string,int> scalar_vars;
    unordered_map<string,int> vector_vars;
    unordered_map<string,int> component_dims;

    void check_not_in_data_mode() ;
    int get_scalar_varid(const string &name) ;
    int get_vector_varid(const string &name) ;
};

/** **/
NetCDFWriter::NetCDFWriter(const string &filename) : ncid(-1), current_time(0)
{
  int ret = nc_create(filename.c_str(), NC_CLOBBER, &ncid);
  if (ret != NC_NOERR)
    flog << "Failed to create NetCDF file: " << nc_strerror(ret);
}

/** **/
NetCDFWriter::~NetCDFWriter() 
{
  if (ncid >= 0) nc_close(ncid);
}

/** **/
void NetCDFWriter::define_dimensions(size_t npoints) 
{
  check_not_in_data_mode();
  nc_def_dim(ncid, "time", NC_UNLIMITED, &dim_time);
  nc_def_dim(ncid, "points", npoints, &dim_points);

  // Coordinate variables
  nc_def_var(ncid, "X", NC_DOUBLE, 1, &dim_points, &var_x);
  nc_def_var(ncid, "Y", NC_DOUBLE, 1, &dim_points, &var_y);
  nc_def_var(ncid, "Z", NC_DOUBLE, 1, &dim_points, &var_z);
}

/** **/
void NetCDFWriter::define_scalar(const string &name) 
{
  check_not_in_data_mode();
  int dims[2] = {dim_time, dim_points};
  int varid;
  nc_def_var(ncid, name.c_str(), NC_DOUBLE, 2, dims, &varid);
  scalar_vars[name] = varid;
}

/** **/
void NetCDFWriter::define_vector(const string &name, int ncomp) 
{
  check_not_in_data_mode();
  int comp_dim;
  auto it = component_dims.find(name);
  if (it == component_dims.end()) {
    nc_def_dim(ncid, (name + "_comp").c_str(), ncomp, &comp_dim);
    component_dims[name] = comp_dim;
  } else {
    comp_dim = it->second;
  }
  int dims[3] = {dim_time, dim_points, comp_dim};
  int varid;
  nc_def_var(ncid, name.c_str(), NC_DOUBLE, 3, dims, &varid);
  vector_vars[name] = varid;
}

/** **/
void NetCDFWriter::define_tensor(const string &name) 
{
  define_vector(name, 9);
}

/** Call this after defining all variables and before any data writing. **/
void NetCDFWriter::end_define() 
{
  int err = nc_enddef(ncid);
  if (err != NC_NOERR)
    flog << "Error ending define mode: " << nc_strerror(err);
  data_mode = true;
}

/** **/
void NetCDFWriter::set_time(size_t t) 
{ current_time = t; }

/** **/
void NetCDFWriter::set_coord(size_t point_idx, double x, double y, double z) 
{
  size_t start[1] = {point_idx}, count[1] = {1};
  nc_put_vara_double(ncid, var_x, start, count, &x);
  nc_put_vara_double(ncid, var_y, start, count, &y);
  nc_put_vara_double(ncid, var_z, start, count, &z);
}

/** **/
void NetCDFWriter::set_coord(size_t point_idx, const Point &p) 
{
  set_coord(point_idx, p(0), p(1), p(2));
}

/** **/
void NetCDFWriter::set_scalar(const string &name, size_t point_idx, double val) 
{
  int varid = get_scalar_varid(name);
  size_t start[2] = {current_time, point_idx}, count[2] = {1,1};
  nc_put_vara_double(ncid, varid, start, count, &val);
}
/** **/
void NetCDFWriter::set_vector(const string &name, size_t point_idx, const vector<double> &vals) 
{
  int varid = get_vector_varid(name);
  size_t start[3] = {current_time, point_idx, 0}, count[3] = {1,1,vals.size()};
  nc_put_vara_double(ncid, varid, start, count, vals.data());
}
/** **/
void NetCDFWriter::set_point(const string &name, size_t point_idx, const vector<double> &vals) 
{
  if (vals.size() != 3) flog << "set_point: expected 3 components";
  set_vector(name, point_idx, vals);
}
/** **/
void NetCDFWriter::set_point(const string &name, size_t point_idx, const Point &p) 
{
  set_vector(name, point_idx, {p(0), p(1), p(2)});
}
/** **/
void NetCDFWriter::set_tensor(const string &name, size_t point_idx, const vector<vector<double>> &vals) 
{
  if (vals.size() != 3 || vals[0].size() != 3 || vals[1].size() != 3 || vals[2].size() != 3) flog << "Expected a 3x3 tensor";

  // Flatten
  vector<double> flat(9);
  for (size_t i=0;i<3;i++)
  for (size_t j=0;j<3;j++)
    flat[i*3+j] = vals[i][j];

  set_vector(name, point_idx, flat);
}

/** **/
void NetCDFWriter::set_tensor(const string &name, size_t point_idx, const RealTensor &T) 
{
  // Point to vector<vector>>
  vector<vector<double>> vals(3, vector<double>(3));
  for (unsigned i=0;i<3;i++)
  for (unsigned j=0;j<3;j++)
    vals[i][j] = T(i,j);

  set_tensor(name, point_idx, vals);
}
/** **/
void NetCDFWriter::define_unit(const string &var_name, const string &units) 
{
  int varid = -1;

  // Determine which map contains the variable
  auto it_scalar = scalar_vars.find(var_name);
  auto it_vector = vector_vars.find(var_name);

  if (it_scalar != scalar_vars.end()) 
  {
    varid = it_scalar->second;
  } else if (it_vector != vector_vars.end()) 
  {
    varid = it_vector->second;
  } else 
  {
    flog << "Variable '" << var_name << "' not defined";
  }

  // Set the "units" attribute
  int ret = nc_put_att_text(ncid, varid, "units", units.size(), units.c_str());
  if (ret != NC_NOERR) 
    flog << "Failed to set units for '" << var_name << "': " << nc_strerror(ret);
}    

/** **/
void NetCDFWriter::check_not_in_data_mode() 
{
  if (data_mode) flog <<  "Cannot define dimensions/variables after writing data";
}
/** **/
int NetCDFWriter::get_scalar_varid(const string &name) 
{
  auto it = scalar_vars.find(name);
  if (it == scalar_vars.end()) flog << "Scalar variable '" << name << "' not defined";
  return it->second;
}
/** **/
int NetCDFWriter::get_vector_varid(const string &name) 
{
  auto it = vector_vars.find(name);
  if (it == vector_vars.end()) flog << "Vector variable '" << name << "' not defined";
  return it->second;
}


/**
 *
 *
 */
int main(int argc, char** argv) {
  HarpyInit init( argc, argv );

  size_t npoints = 4;
  NetCDFWriter writer("points_libmesh.nc");

  // 1️⃣ Define dimensions & variables first
  writer.define_dimensions(npoints);
  writer.define_scalar("Temperature");
  writer.define_unit( "Temperature", "C" );
  writer.define_vector("Velocity", 3);
  writer.define_tensor("Stress");
  writer.end_define();

  // 2️⃣ Set coordinates
  vector<Point> points = {
    Point(0.0, 0.0, 0.0),
    Point(1.0, 0.0, 0.0),
    Point(0.0, 1.0, 0.0),
    Point(1.0, 1.0, 1.0)
  };

  // 3️⃣ Write timestep 0
  writer.set_time(0);
  for (size_t i = 0; i < npoints; ++i) 
  {
    writer.set_coord(i, points[i]);
    writer.set_scalar("Temperature", i, 300.0+i);
    writer.set_point("Velocity", i, points[i]); // Velocity = point coordinates
    RealTensor T;
    T(0,0)=1.0; T(0,1)=0.1; T(0,2)=0.2;
    T(1,0)=0.1; T(1,1)=1.1; T(1,2)=0.3;
    T(2,0)=0.2; T(2,1)=0.3; T(2,2)=1.2;
    writer.set_tensor("Stress", i, T);
  }

  cout << "NetCDF file 'points_libmesh.nc' written successfully!" << endl;


  return 0;
}

