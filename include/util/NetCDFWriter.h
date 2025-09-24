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
enum class NC_TYPE   { SCALAR, VEC3, TEN9 };
enum class NC_PARAM  { 
  TEMPERATURE, DELTA_T,
  PRESSURE, DELTA_P,
  VELOCITY, 
  DENSITY, VISCOSITY, 
  S1, S1_MAG ,
  S2, S2_MAG ,
  S3, S3_MAG ,
  SIGTOT,
  VP_STRAIN ,
  VP_STRAIN_RATE ,
  INVAR_P_EFF, INVAR_Q,
};

struct DataInfo { 
  NC_TYPE type;
  string name;
  string unit;
  string description;
  int vid;
};

extern map<NC_PARAM, DataInfo> PARAMS;

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

    map<NC_PARAM, vector<float>> buffer;

  public:    

    uint n_points;

    // File state
    string filename;
    bool file_created;
    bool define_mode;

    uint curr_pt, curr_time;
    void set_curr_pt(uint i)   { curr_pt = i; }
    void set_curr_time(uint i) { curr_time = i; }
    void add_timestep(uint i, double time);

    void flush( bool sync=true );

    NetCDFWriter( string filename_ );
    ~NetCDFWriter();

    // File management
    void init( uint n_points_ );
    void close_file();

    // Variable definition methods
    void add( NC_PARAM var );
    void finish_definitions();

    // Coordinate and time setting methods
    void set_coords(const vector<Point>& points);
    void set_time_collective(const vector<double>& times);

    // Data setting methods
    void set_value(NC_PARAM var, double value);
    void set_value(NC_PARAM var, const Point& vector);
    void set_value(NC_PARAM var, const RealTensorValue& tensor);

    // Zero out the buffers
    void reset_buffer() { for ( auto & [ p, vec ] : buffer ) vec.assign(vec.size(), 0.0f); }

  private:
    inline DataInfo & get_di( NC_PARAM var );
    void setup_coordinate_variables();
    void add_global_attributes();
    void set_collective_access(int varid);
    void set_independent_access(int varid);
};

std::ostream& operator<<(std::ostream& os, const NC_TYPE& type);
std::ostream& operator<<(std::ostream& os, const NC_PARAM& param);


/** **/
inline DataInfo & NetCDFWriter::get_di( NC_PARAM var )
{
  if ( ! PARAMS.count( var ) ) flog << "Unknown variable " << var << ".";
  return PARAMS.at(var);
}

} // ns
