
#include "util/NetCDFWriter.h"

namespace util {

  /**
   *   Debug stuff
   */
void debug_check_time_dimension(int ncid, int time_dim_id, int rank) {
    size_t time_len;
    int retval;
    
    // Sync first to ensure we see latest changes
    nc_sync(ncid);
    
    // Query the dimension length
    if ((retval = nc_inq_dimlen(ncid, time_dim_id, &time_len))) {
        printf("[Rank %d] Error getting time dimension: %s\n", 
               rank, nc_strerror(retval));
    } else {
        printf("[Rank %d] Time dimension length: %zu\n", rank, time_len);
    }
    
    // Flush output immediately
    fflush(stdout);
}

map<NC_PARAM, DataInfo> PARAMS = {
    { NC_PARAM::TEMPERATURE,      {  NC_TYPE::SCALAR, "Temperature"       , "K",    "", -1} },
    { NC_PARAM::PRESSURE,         {  NC_TYPE::SCALAR, "Pressure"          , "Pa",   "", -1} },
    { NC_PARAM::DELTA_T    ,      {  NC_TYPE::SCALAR, "Delta_T"           , "K",    "", -1} },
    { NC_PARAM::DELTA_P,          {  NC_TYPE::SCALAR, "Delta_P"           , "Pa",   "", -1} },
    { NC_PARAM::S1,               {  NC_TYPE::VEC3,   "S1"                , "Pa",   "", -1} },
    { NC_PARAM::S1_MAG,           {  NC_TYPE::SCALAR, "S1 Magnitude"      , "Pa",   "", -1} },
    { NC_PARAM::S2,               {  NC_TYPE::VEC3,   "S2"                , "Pa",   "", -1} },
    { NC_PARAM::S2_MAG,           {  NC_TYPE::SCALAR, "S2 Magnitude"      , "Pa",   "", -1} },
    { NC_PARAM::S3,               {  NC_TYPE::VEC3,   "S3"                , "Pa",   "", -1} },
    { NC_PARAM::S3_MAG,           {  NC_TYPE::SCALAR, "S3 Magnitude"      , "Pa",   "", -1} },
    { NC_PARAM::INVAR_P_EFF,      {  NC_TYPE::SCALAR, "Invariant P(eff)"  , "Pa",   "", -1} },
    { NC_PARAM::INVAR_Q,          {  NC_TYPE::SCALAR, "Invariant Q"       , "Pa",   "", -1} },
    { NC_PARAM::VELOCITY,         {  NC_TYPE::VEC3,   "Velocity"          , "m/s",  "", -1} },
    { NC_PARAM::SIGTOT,           {  NC_TYPE::TEN9,   "Total Stress"      , "Pa",   "", -1} },
    { NC_PARAM::VP_STRAIN,        {  NC_TYPE::TEN9,   "VP Strain"         , "",     "", -1} },
    { NC_PARAM::VP_STRAIN_RATE,   {  NC_TYPE::TEN9,   "VP Strain Rate"    , "1/s",  "", -1} }
};

/** **/
NetCDFWriter::NetCDFWriter( string filename_ )
  : world(), n_points(0), filename(filename_), 
  file_created(false), define_mode(false),
  vec3_dimid(0), ten9_dimid(0), curr_time(0), curr_pt(0) { }

  /** **/
  NetCDFWriter::~NetCDFWriter() {
    if (file_created) close_file();
  }

/** **/
void NetCDFWriter::init( uint n_points_ ) 
{
  if (file_created) flog << "File already created";

  n_points = n_points_;
  MPI_Info info = MPI_INFO_NULL;

  // Create parallel NetCDF file
  CHECK_NC( 
      nc_create_par(filename.c_str(), NC_NETCDF4 | NC_MPIIO, world, info, &ncid)
  );

  // Define dimensions
  CHECK_NC( nc_def_dim(ncid, "time", NC_UNLIMITED, &time_dimid)   );
  CHECK_NC( nc_def_dim(ncid, "point_idx", n_points, &point_dimid) );
  CHECK_NC( nc_def_dim(ncid, "vec3_comp", 3, &vec3_dimid)         );

  setup_coordinate_variables();
  add_global_attributes();

  file_created = true;
  define_mode = true;

  if (world.rank() == 0) {
    ilog << "Created NetCDF file: " << filename;
    ilog << "Dimensions: time(UNLIMITED), point_idx(" << n_points << ")";
  }
}

/** **/
void NetCDFWriter::setup_coordinate_variables() 
{
  // Define coordinate variables
  CHECK_NC(nc_def_var(ncid, "time", NC_FLOAT, 1, &time_dimid, &time_varid));

  int dimids[2] = {point_dimid, vec3_dimid};
  CHECK_NC(nc_def_var(ncid, "Coord", NC_FLOAT, 2, dimids, &coord_varid));
  string comp = "x y z";
  CHECK_NC(nc_put_att_text(ncid, coord_varid, "components", comp.size(), comp.c_str()));    

  // Add coordinate attributes
  CHECK_NC(nc_put_att_text(ncid, time_varid, "units", 1, "s"));
  CHECK_NC(nc_put_att_text(ncid, coord_varid, "units", 1, "m"));

  // Collective access to coord and time (all process write in sync)
  set_collective_access( coord_varid );
  set_collective_access( time_varid );
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
void NetCDFWriter::add( NC_PARAM var )
{
  if (!define_mode) flog << "Not in define mode";
  auto & di = get_di(var);

  vector<int> dimids;
  if ( di.type == NC_TYPE::SCALAR ) dimids = {time_dimid, point_dimid};
  if ( di.type == NC_TYPE::VEC3 ) dimids = {time_dimid, point_dimid, vec3_dimid};
  if ( di.type == NC_TYPE::TEN9 )   
  {
    // Create the dimension, if not existing
    if ( ! ten9_dimid )  CHECK_NC(nc_def_dim(ncid, "ten9_comp", 9, &ten9_dimid));

    dimids = {time_dimid, point_dimid, ten9_dimid};
  }

  CHECK_NC(nc_def_var(ncid, di.name.c_str(), NC_FLOAT, dimids.size(), dimids.data(), &di.vid));

  CHECK_NC(nc_put_att_text(ncid, di.vid, "units", di.unit.length(), di.unit.c_str()));

  // Useful annotations
  if ( di.type == NC_TYPE::TEN9 ) 
  {
    string comp = "xx xy xz yx yy yz zx zy zz";
    CHECK_NC(nc_put_att_text(ncid, di.vid, "components", comp.size(), comp.c_str()));    
  }

  if ( di.type == NC_TYPE::VEC3 ) 
  {
    string comp = "x y z";
    CHECK_NC(nc_put_att_text(ncid, di.vid, "components", comp.size(), comp.c_str()));    
  }

  // Independent access: all processors write independently (not sync'ed)
  set_independent_access( di.vid );

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
void NetCDFWriter::set_coords(const vector<Point>& points) 
{
  if (define_mode) flog << "Still in define mode - call finish_definitions() first";

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

/**
 *      This is a collective task
 */
void NetCDFWriter::add_timestep(uint i, double time)
{
  if (define_mode) flog << "Still in define mode - call finish_definitions() first";

  // Write time data collectively (all processes write same data)
  curr_time = i;
  luint start[1] = {curr_time};
  luint count[1] = {1};

  // Extend the dimesion of the time variable
  set_collective_access( time_varid );
  CHECK_NC(nc_put_vara_double(ncid, time_varid, start, count, &time));

  // Extend the dimension of all variables (needs to be collective!)
  for ( auto & [ p, di ] : PARAMS )
  if ( di.vid >= 0 )
  {
    set_collective_access( di.vid );
    if ( di.type == NC_TYPE::SCALAR ) set_value( p, 0 );
    if ( di.type == NC_TYPE::VEC3 )   set_value( p, Point() );
    if ( di.type == NC_TYPE::TEN9 )   set_value( p, RealTensorValue() );
    set_independent_access( di.vid );
  }

  // Ensures all processes are synchronized  
  world.barrier(); 
}

///** **/
//void NetCDFWriter::set_time_collective(const vector<double>& times) 
//{
//  if (define_mode) flog << "Still in define mode - call finish_definitions() first";

//  if (times.size() != n_times) flog << "Times vector size doesn't match n_times";

//  // Set collective access for time variable
//  set_collective_access(time_varid);

//  // Write time data collectively (all processes write same data)
//  luint start[1] = {0};
//  luint count[1] = {n_times};

//  CHECK_NC(nc_put_vara_double(ncid, time_varid, start, count, times.data()));

//  if (world.rank() == 0) ilog << "Set time collectively";
//}

/** **/
void NetCDFWriter::set_value(NC_PARAM var, double value) 
{
  if (define_mode) flog << "Still in define mode - call finish_definitions() first";
  auto & di = get_di(var);

  luint start[2] = { curr_time, curr_pt };
  luint count[2] = { 1, 1 };

  CHECK_NC(nc_put_vara_double(ncid, di.vid, start, count, &value));
}

/** **/
void NetCDFWriter::set_value(NC_PARAM var, const Point& p) 
{
  if (define_mode) flog << "Still in define mode - call finish_definitions() first";
  auto & di = get_di(var);

  luint start[3] = { curr_time, curr_pt, 0};
  luint count[3] = {1, 1, 3};

  std::vector<double> vec = {p(0), p(1), p(2)};
  CHECK_NC( nc_put_vara_double( ncid, di.vid, start, count, vec.data() ) );
}

/** **/
void NetCDFWriter::set_value(NC_PARAM var, const RealTensorValue& tensor) 
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

/** **/
ostream& operator<<(ostream& os, const NC_TYPE& type) {
    static const unordered_map<NC_TYPE, string> type_map = {
        {NC_TYPE::SCALAR, "SCALAR"},
        {NC_TYPE::VEC3,   "VEC3"},
        {NC_TYPE::TEN9,   "TEN9"}
    };
    
    auto it = type_map.find(type);
    if (it != type_map.end()) return os << it->second;
    return os << "Unknown NC_TYPE(" << static_cast<int>(type) << ")";
}
/** **/
ostream& operator<<(ostream& os, const NC_PARAM& param) {
    static const unordered_map<NC_PARAM, string> param_map = {
        {NC_PARAM::DELTA_T,           "DELTA_T"},
        {NC_PARAM::DELTA_P,           "DELTA_P"},
        {NC_PARAM::INVAR_P_EFF,       "INVAR_P_EFF"},
        {NC_PARAM::INVAR_Q,           "INVAR_Q"},
        {NC_PARAM::TEMPERATURE,       "TEMPERATURE"},
        {NC_PARAM::PRESSURE,          "PRESSURE"},
        {NC_PARAM::S1_MAG,            "S1_MAG"},
        {NC_PARAM::S1,                "S1"},
        {NC_PARAM::S2_MAG,            "S2_MAG"},
        {NC_PARAM::S2,                "S2"},
        {NC_PARAM::S3_MAG,            "S3_MAG"},
        {NC_PARAM::S3,                "S3"},
        {NC_PARAM::VELOCITY,          "VELOCITY"},
        {NC_PARAM::DENSITY,           "DENSITY"},
        {NC_PARAM::VISCOSITY,         "VISCOSITY"},
        {NC_PARAM::SIGTOT,            "SIGTOT"},
        {NC_PARAM::VP_STRAIN,         "VP_STRAIN"},
        {NC_PARAM::VP_STRAIN_RATE,    "VP_STRAIN_RATE"}
    };
    
    auto it = param_map.find(param);
    if (it != param_map.end()) return os << it->second;
    return os << "Unknown NC_PARAM(" << static_cast<int>(param) << ")";
}

} // ns
