#pragma once

#include "harpy/Global.h"

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

namespace util {

using libMesh::Point;
using libMesh::RealTensorValue;

namespace mpi = boost::mpi;


// Variable types enum
enum VAR_TYPE { SCALAR, VEC3, TEN9 };
enum NC_VAR   { TEMPERATURE, PRESSURE, VELOCITY, STRESS, DENSITY, VISCOSITY };
struct DataInfo { 
  VAR_TYPE type;
  string name;
  string unit;
  string description;
  int vid;
};

extern map<NC_VAR, DataInfo> PARAMS;

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
    int coord_varid, time_varid, ten9_varid;

  public:    

    uint n_times;
    uint n_points;

    // File state
    bool file_created;
    bool define_mode;

    uint curr_pt, curr_time;
    void set_curr_pt(uint i) { curr_pt = i; }
    void set_curr_time(uint i) { curr_time = i; }

    explicit NetCDFWriter(uint n_times, uint n_points);
    ~NetCDFWriter();

    // File management
    void create_file(const string& filename);
    void close_file();

    // Variable definition methods
    void define_variable( NC_VAR var );
    void finish_definitions();

    // Coordinate and time setting methods
    void set_coordinates_collective(const vector<Point>& points);
    void set_time_collective(const vector<double>& times);

    // Data setting methods
    void set_value(NC_VAR var, double value);
    void set_value(NC_VAR var, const Point& vector);
    void set_value(NC_VAR var, const RealTensorValue& tensor);

  private:
    inline DataInfo & get_di( NC_VAR var );
    void setup_coordinate_variables();
    void add_global_attributes();
    void set_collective_access(int varid);
    void set_independent_access(int varid);
};


/** **/
inline DataInfo & NetCDFWriter::get_di( NC_VAR var )
{
  if ( ! PARAMS.count( var ) ) flog << "Unknown variable " << var << ".";
  return PARAMS.at(var);
}

} // ns
