
#include "util/NetCDFWriter.h"
#include "util/OutputOperators.h"

namespace mpi = boost::mpi;

namespace util {


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
    { NC_PARAM::VP_STRAIN_RATE,   {  NC_TYPE::TEN9,   "VP Strain Rate"    , "1/s",  "", -1} },
    { NC_PARAM::U,                {  NC_TYPE::VEC3,   "U"                 , "m",    "", -1} },
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
//      nc_create_par(filename.c_str(), NC_NETCDF4 | NC_MPIIO, world, info, &ncid)
      nc_create_par(filename.c_str(), NC_NETCDF4 , world, info, &ncid)
  );

  // Increase chunk size to get more performance in NFS systems
  nc_set_chunk_cache(512*1024*1024, 4133, 0.75);
  {
  size_t size, nelems;
  float preemption;
  nc_get_chunk_cache(&size, &nelems, &preemption);
//  ilog << "Cache: " << size << " bytes, " << nelems << " slots, " << preemption << " preemption";
  }

  // Define dimensions
  CHECK_NC( nc_def_dim(ncid, "time", NC_UNLIMITED, &time_dimid)   );
  CHECK_NC( nc_def_dim(ncid, "point_idx", n_points, &point_dimid) );
  CHECK_NC( nc_def_dim(ncid, "vec3_comp", 3, &vec3_dimid)         );

  setup_coordinate_variables();
  add_global_attributes();

  file_created = true;
  define_mode = true;
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

  uint n_floats=1;
  vector<int> dimids;
  if ( di.type == NC_TYPE::SCALAR ) 
  {
    dimids = {time_dimid, point_dimid};
    n_floats = 1;
  }

  if ( di.type == NC_TYPE::VEC3 ) 
  {
    n_floats = 3;
    dimids = {time_dimid, point_dimid, vec3_dimid};
  }

  if ( di.type == NC_TYPE::TEN9 )   
  {
    n_floats = 9;

    // Create the dimension, if not existing
    if ( ! ten9_dimid )  CHECK_NC(nc_def_dim(ncid, "ten9_comp", 9, &ten9_dimid));

    dimids = {time_dimid, point_dimid, ten9_dimid};
  }

  // Allocate buffer
  buffer[var] = std::vector<float>();
  buffer[var].resize( n_points * n_floats );

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

  // Variables are independent (only assigned by rank=0)
  set_independent_access( di.vid );
}

/** **/
void NetCDFWriter::finish_definitions() 
{
  if (!define_mode) flog << "Not in define mode";
  CHECK_NC(nc_enddef(ncid));
  define_mode = false;
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

  reset_buffer();

  // Assign an empty time entry in collective mode to ensure the dimension is properly extended
  for (auto& [var, data] : buffer) 
  {
    DataInfo & di = get_di( var );
    vector<luint> start, count;
    if ( di.type == NC_TYPE::SCALAR ) 
    {
      start = { curr_time  , 0          };
      count = { 1          , n_points    };
      ASSERT( data.size() == n_points , "[scalar] Data size not matching. " << data.size() << " != " << n_points << "." );
    } 
    else if ( di.type == NC_TYPE::VEC3 )   
    {
      start = { curr_time  , 0        , 0 };
      count = { 1          , n_points  , 3 };
      ASSERT( data.size() == n_points*3 , "[vec3] Data size not matching. " << data.size() << " != " << n_points << "." );
    } 
    else if ( di.type == NC_TYPE::TEN9 )
    {
      start = { curr_time , 0       , 0 };
      count = { 1         , n_points , 9 };
      ASSERT( data.size() == n_points*9 , "[ten9] Data size not matching. " << data.size() << " != " << n_points << "." );
    } else flog << "Should not be here...";

    //
    set_collective_access( di.vid ); // All vars are always collective
    CHECK_NC( nc_put_vara_float(ncid, di.vid, start.data(), count.data(), data.data() ));
    set_independent_access( di.vid );
    //
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

  ASSERT( curr_pt < buffer[var].size(), "Trying to write outside the buffer. (" << curr_pt << " >= " << buffer[var].size() << " ) for var " << var << "" );
  buffer[var][curr_pt] = value;
//  CHECK_NC(nc_put_vara_double(ncid, di.vid, start, count, &value));
}

/** **/
void NetCDFWriter::set_value(NC_PARAM var, const Point& p) 
{
  if (define_mode) flog << "Still in define mode - call finish_definitions() first";
  auto & di = get_di(var);

  luint start[3] = { curr_time, curr_pt, 0};
  luint count[3] = {1, 1, 3};

  std::vector<double> vec = {p(0), p(1), p(2)};

  ASSERT( 3*curr_pt+2 < buffer[var].size(), "Trying to write outside the buffer. (" << 3*curr_pt+2 << " >= " << buffer[var].size() << " ) for var " << var << "" );
  for ( uint i = 0; i < 3 ; i++ ) buffer[var][3*curr_pt+i] = p(i);
//  CHECK_NC( nc_put_vara_double( ncid, di.vid, start, count, vec.data() ) );
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

  ASSERT( 9*curr_pt+2 < buffer[var].size(), "Trying to write outside the buffer. (" << 9*curr_pt+2 << " >= " << buffer[var].size() << " ) for var " << var << "" );
  for ( uint i = 0; i < 9 ; i++ ) buffer[var][9*curr_pt+i] = flat[i];
//  CHECK_NC( nc_put_vara_double( ncid, di.vid, start, count, flat.data() ) );
}

/** Flushes the buffer into the netcdf **/
void NetCDFWriter::flush( bool sync )
{
  for (auto& [var, data] : buffer) 
  {
    ASSERT( !data.empty() , "Data is empty for var " << var );

    // Perform element-wise sum reduction
    std::vector<float> result;
    if (world.rank() == 0) result.resize(data.size());
    MPI_Reduce(data.data(),           // send buffer
        world.rank() == 0 ? result.data() : nullptr,  // receive buffer (only on rank 0)
        data.size(),           // number of elements
        MPI_FLOAT,            // data type
        MPI_SUM,              // reduction operation
        0,                    // root process (rank 0)
        MPI_COMM_WORLD);      // communicator

    DataInfo & di = get_di( var );

//    checkIfUnlimited( ncid, time_dimid , "Time" );
//    checkVariableDimensions( ncid, di.vid, start, count );

    // Write only at rank=0 (variables must be INDEPENDENT)
    if ( world.rank() == 0 )
    {
      vector<luint> start, count;
      if ( di.type == NC_TYPE::SCALAR ) 
      {
        start = { curr_time  , 0          };
        count = { 1          , n_points    };
        ASSERT( result.size() == n_points , "[scalar] Data size not matching. " << result.size() << " != " << n_points << "." );
      } 
      else if ( di.type == NC_TYPE::VEC3 )   
      {
        start = { curr_time  , 0        , 0 };
        count = { 1          , n_points  , 3 };
        ASSERT( result.size() == n_points*3 , "[vec3] Data size not matching. " << result.size() << " != " << n_points << "." );
      } 
      else if ( di.type == NC_TYPE::TEN9 )
      {
        start = { curr_time , 0       , 0 };
        count = { 1         , n_points , 9 };
        ASSERT( result.size() == n_points*9 , "[ten9] Data size not matching. " << result.size() << " != " << n_points << "." );
      } else flog << "Should not be here...";

      CHECK_NC( nc_put_vara_float(ncid, di.vid, start.data(), count.data(), result.data() ));
    }
  } 

  reset_buffer();
  if ( sync ) nc_sync(ncid); 
}

/** **/
void NetCDFWriter::set_collective_access(int varid) { CHECK_NC(nc_var_par_access(ncid, varid, NC_COLLECTIVE)); }
void NetCDFWriter::set_independent_access(int varid) { CHECK_NC(nc_var_par_access(ncid, varid, NC_INDEPENDENT)); }
void NetCDFWriter::close_file() 
{
  if (! file_created) return;
  CHECK_NC(nc_close(ncid));
  file_created = false;
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
        {NC_PARAM::VP_STRAIN_RATE,    "VP_STRAIN_RATE"},
        {NC_PARAM::U,                 "U"},
    };
    
    auto it = param_map.find(param);
    if (it != param_map.end()) return os << it->second;
    return os << "Unknown NC_PARAM(" << static_cast<int>(param) << ")";
}

} // ns

//  /**
//   *   Debug stuff
//   */
//void debug_check_time_dimension(int ncid, int time_dim_id, int rank) {
//    size_t time_len;
//    int retval;
//    
//    // Sync first to ensure we see latest changes
//    nc_sync(ncid);
//    
//    // Query the dimension length
//    if ((retval = nc_inq_dimlen(ncid, time_dim_id, &time_len))) {
//        printf("[Rank %d] Error getting time dimension: %s\n", 
//               rank, nc_strerror(retval));
//    } else {
//        printf("[Rank %d] Time dimension length: %zu\n", rank, time_len);
//    }
//    
//    // Flush output immediately
//    fflush(stdout);
//}
//void checkIfUnlimited(int ncid, int dimid, const char* dimname) {
//    int unlimdimid;
//    nc_inq_unlimdim(ncid, &unlimdimid);

//    size_t dimsize;
//    nc_inq_dimlen(ncid, dimid, &dimsize);

//    std::cout << "Dimension '" << dimname << "' (id=" << dimid << "):" << std::endl;
//    std::cout << "  Current size: " << dimsize << std::endl;
//    std::cout << "  Unlimited dim ID: " << unlimdimid << std::endl;
//    std::cout << "  Is unlimited: " << (dimid == unlimdimid ? "YES" : "NO") << std::endl;
//}

//bool checkVariableDimensions(int ncid, int varid, const std::vector<luint>& start, const std::vector<luint>& count) {
//    int ndims;
//    int dimids[NC_MAX_VAR_DIMS];

//    // Get unlimited dimension ID
//    int unlimdimid;
//    nc_inq_unlimdim(ncid, &unlimdimid);

//    char varname[NC_MAX_NAME+1];
//    nc_inq_var(ncid, varid, varname, nullptr, &ndims, dimids, nullptr);

//    ilog1 << "Variable " << varname << " has " << ndims << " dimensions:" << std::endl;

//    bool needs_collective = false;
//    for (int i = 0; i < ndims; i++) {
//        char dimname[NC_MAX_NAME+1];
//        size_t current_size;

//        nc_inq_dim(ncid, dimids[i], dimname, &current_size);

//        // CORRECT WAY: Compare dimension IDs, not sizes
//        bool is_unlimited = (dimids[i] == unlimdimid);

//        size_t needed_size = start[i] + count[i];
//        cout << "i:" << i << " start[i]:" << start[i] << " count[i]:" << count[i] << endl;

//        ilog1 << "  Dim[" << i << "] '" << dimname << "': current=" << current_size
//                  << " needed=" << needed_size << " unlimited=" << is_unlimited << std::endl;

//        if (is_unlimited) {
//            ilog1 << "    INFO: Unlimited dimension - requires COLLECTIVE mode" << std::endl;
//            needs_collective = true;
//        } else if (needed_size > current_size) {
//            ilog1 << "    ERROR: Would exceed fixed dimension!" << std::endl;
//        }
//    }

//    if (needs_collective) {
//        ilog1 << "Variable requires COLLECTIVE mode (has unlimited dimension)" << std::endl;
//    } else {
//        ilog1 << "Variable can use INDEPENDENT mode (no unlimited dimensions)" << std::endl;
//    }

//    return needs_collective;
//}

