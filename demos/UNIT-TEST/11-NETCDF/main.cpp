#include <netcdf.h>
#include <iostream>
#include <vector>
#include <string>
#include <unordered_map>
#include <stdexcept>

using namespace std;

class NetCDFWriter {
public:
    NetCDFWriter(const string &filename)
        : ncid(-1), current_time(0), current_x(0), current_y(0), current_z(0)
    {
        if (nc_create(filename.c_str(), NC_CLOBBER, &ncid))
            throw runtime_error("Error creating NetCDF file");
    }

    ~NetCDFWriter() {
        if (ncid >= 0) nc_close(ncid);
    }

    // Define dimensions: time unlimited, x/y/z fixed
    void define_dimensions(int nx, int ny, int nz) {
        nc_def_dim(ncid, "time", NC_UNLIMITED, &dim_time);
        nc_def_dim(ncid, "x", nx, &dim_x);
        nc_def_dim(ncid, "y", ny, &dim_y);
        nc_def_dim(ncid, "z", nz, &dim_z);
    }

    void define_scalar(const string &name) {
        int dims[4] = {dim_time, dim_x, dim_y, dim_z};
        int varid;
        nc_def_var(ncid, name.c_str(), NC_FLOAT, 4, dims, &varid);
        scalar_vars[name] = varid;
    }

    void define_vector(const string &name, int components) {
        int comp_dim;
        nc_def_dim(ncid, (name + "_comp").c_str(), components, &comp_dim);
        int dims[5] = {dim_time, dim_x, dim_y, dim_z, comp_dim};
        int varid;
        nc_def_var(ncid, name.c_str(), NC_FLOAT, 5, dims, &varid);
        vector_vars[name] = varid;
    }

    void define_tensor(const string &name, int components) {
        define_vector(name, components); // treat as vector with components
    }

    void end_define() {
        if (nc_enddef(ncid))
            throw runtime_error("Error ending define mode");
    }

    // Set current indices
    void set_time(int t) { current_time = t; }
    void set_point(int x, int y, int z) { current_x = x; current_y = y; current_z = z; }

    // Write scalar
    void set_var(const string &name, float val) {
        auto varid = scalar_vars.at(name);
        size_t start[4] = {static_cast<size_t>(current_time),
                           static_cast<size_t>(current_x),
                           static_cast<size_t>(current_y),
                           static_cast<size_t>(current_z)};
        size_t count[4] = {1,1,1,1};
        nc_put_vara_float(ncid, varid, start, count, &val);
    }

    // Write vector/tensor
    void set_var(const string &name, const vector<float> &vals) {
        auto varid = vector_vars.at(name);
        size_t start[5] = {static_cast<size_t>(current_time),
                           static_cast<size_t>(current_x),
                           static_cast<size_t>(current_y),
                           static_cast<size_t>(current_z),
                           0};
        size_t count[5] = {1,1,1,1, vals.size()};
        nc_put_vara_float(ncid, varid, start, count, vals.data());
    }

    // Write a single 3x3 tensor at the current point
    void set_var(const std::string &name, const std::vector<std::vector<float>> &vals) {
      // Check size
      if (vals.size() != 3 || vals[0].size() != 3 ||
          vals[1].size() != 3 || vals[2].size() != 3)
        throw std::runtime_error("set_var: expected a 3x3 tensor");

      // Flatten 3x3 tensor into 1D vector
      std::vector<float> flat(9);
      for (size_t i = 0; i < 3; ++i)
        for (size_t j = 0; j < 3; ++j)
          flat[i*3 + j] = vals[i][j];

      // Reuse the existing vector version
      set_var(name, flat);
    }

private:
    int ncid;
    int dim_time, dim_x, dim_y, dim_z;
    int current_time, current_x, current_y, current_z;

    unordered_map<string,int> scalar_vars;
    unordered_map<string,int> vector_vars;
};

// Example usage
int main() {
  int nx = 3, ny = 4, nz = 2;
  int nt = 2;

  NetCDFWriter obj("data.nc");

  obj.define_dimensions(nx, ny, nz);
  obj.define_scalar("T");
  obj.define_tensor("sigma", 9); // 9-component tensor
  obj.define_vector("vel", 3);   // 3-component vector
  obj.end_define();

  // Loop over time and points
  for (int t=0; t<nt; ++t) {
    obj.set_time(t);
    for (int x=0; x<nx; ++x)
      for (int y=0; y<ny; ++y)
        for (int z=0; z<nz; ++z)
        {
          obj.set_point(x,y,z);

          obj.set_var("T", 100.0f + t);

          vector<vector<float>> sigma(3, vector<float>(3));
          for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
              sigma[i][j] = static_cast<float>(i * j + z * t);

          // Write to NetCDF using the wrapper
          obj.set_var("sigma", sigma);

          vector<float> vel{0.1f*t, 0.2f*t, 0.3f*t};
          obj.set_var("vel", vel);
        }
  }

  cout << "NetCDF file with unlimited time and fixed x/y/z created!" << endl;
}

