#include <netcdf.h>
#include <vector>
#include <string>
#include <unordered_map>
#include <stdexcept>
#include <iostream>
#include "libmesh/point.h"
#include "libmesh/tensor_value.h"

using namespace std;
using namespace libMesh;

class NetCDFWriter {
public:
    NetCDFWriter(const string &filename)
        : ncid(-1), current_time(0)
    {
        int ret = nc_create(filename.c_str(), NC_CLOBBER, &ncid);
        if (ret != NC_NOERR)
            throw runtime_error("Failed to create NetCDF file: " + string(nc_strerror(ret)));
    }

    ~NetCDFWriter() {
        if (ncid >= 0) nc_close(ncid);
    }

    // Define dimensions
    void define_dimensions(size_t npoints) {
        check_not_in_data_mode();
        nc_def_dim(ncid, "time", NC_UNLIMITED, &dim_time);
        nc_def_dim(ncid, "points", npoints, &dim_points);

        // Coordinate variables
        nc_def_var(ncid, "X", NC_DOUBLE, 1, &dim_points, &var_x);
        nc_def_var(ncid, "Y", NC_DOUBLE, 1, &dim_points, &var_y);
        nc_def_var(ncid, "Z", NC_DOUBLE, 1, &dim_points, &var_z);
    }

    // Define scalar [time, points]
    void define_scalar(const string &name) {
        check_not_in_data_mode();
        int dims[2] = {dim_time, dim_points};
        int varid;
        nc_def_var(ncid, name.c_str(), NC_DOUBLE, 2, dims, &varid);
        scalar_vars[name] = varid;
    }

    // Define vector [time, points, components]
    void define_vector(const string &name, int ncomp) {
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

    // Define tensor as 3x3 vector
    void define_tensor(const string &name) {
        define_vector(name, 9);
    }

    void end_define() {
        int err = nc_enddef(ncid);
        if (err != NC_NOERR)
            throw runtime_error(string("Error ending define mode: ") + nc_strerror(err));
        data_mode = true;
    }

    void set_time(size_t t) { current_time = t; }

    // Write coordinates
    void set_coord(size_t point_idx, double x, double y, double z) {
        size_t start[1] = {point_idx}, count[1] = {1};
        nc_put_vara_double(ncid, var_x, start, count, &x);
        nc_put_vara_double(ncid, var_y, start, count, &y);
        nc_put_vara_double(ncid, var_z, start, count, &z);
    }

    void set_coord(size_t point_idx, const Point &p) {
        set_coord(point_idx, p(0), p(1), p(2));
    }

    // Scalars
    void set_scalar(const string &name, size_t point_idx, double val) {
        int varid = get_scalar_varid(name);
        size_t start[2] = {current_time, point_idx}, count[2] = {1,1};
        nc_put_vara_double(ncid, varid, start, count, &val);
    }

    // Vectors
    void set_vector(const string &name, size_t point_idx, const vector<double> &vals) {
        int varid = get_vector_varid(name);
        size_t start[3] = {current_time, point_idx, 0}, count[3] = {1,1,vals.size()};
        nc_put_vara_double(ncid, varid, start, count, vals.data());
    }

    void set_point(const string &name, size_t point_idx, const vector<double> &vals) {
        if (vals.size() != 3) throw runtime_error("set_point: expected 3 components");
        set_vector(name, point_idx, vals);
    }

    void set_point(const string &name, size_t point_idx, const Point &p) {
        set_vector(name, point_idx, {p(0), p(1), p(2)});
    }

    // Tensors
    void set_tensor(const string &name, size_t point_idx, const vector<vector<double>> &vals) {
        if (vals.size() != 3 || vals[0].size() != 3 || vals[1].size() != 3 || vals[2].size() != 3)
            throw runtime_error("Expected a 3x3 tensor");
        vector<double> flat(9);
        for (size_t i=0;i<3;i++)
            for (size_t j=0;j<3;j++)
                flat[i*3+j] = vals[i][j];
        set_vector(name, point_idx, flat);
    }

    void set_tensor(const string &name, size_t point_idx, const RealTensor &T) {
        vector<vector<double>> vals(3, vector<double>(3));
        for (unsigned i=0;i<3;i++)
            for (unsigned j=0;j<3;j++)
                vals[i][j] = T(i,j);
        set_tensor(name, point_idx, vals);
    }

    void define_unit(const string &var_name, const string &units) {
        int varid = -1;

        // Determine which map contains the variable
        auto it_scalar = scalar_vars.find(var_name);
        auto it_vector = vector_vars.find(var_name);

        if (it_scalar != scalar_vars.end()) {
            varid = it_scalar->second;
        } else if (it_vector != vector_vars.end()) {
            varid = it_vector->second;
        } else {
            throw runtime_error("Variable '" + var_name + "' not defined");
        }

        // Set the "units" attribute
        int ret = nc_put_att_text(ncid, varid, "units", units.size(), units.c_str());
        if (ret != NC_NOERR) {
            throw runtime_error("Failed to set units for '" + var_name + "': " + nc_strerror(ret));
        }
    }    

private:
    int ncid;
    int dim_time, dim_points;
    int var_x, var_y, var_z;
    size_t current_time;
    bool data_mode = false;

    unordered_map<string,int> scalar_vars;
    unordered_map<string,int> vector_vars;
    unordered_map<string,int> component_dims;

    void check_not_in_data_mode() {
        if (data_mode) throw runtime_error("Cannot define dimensions/variables after writing data");
    }

    int get_scalar_varid(const string &name) {
        auto it = scalar_vars.find(name);
        if (it == scalar_vars.end()) throw runtime_error("Scalar variable '" + name + "' not defined");
        return it->second;
    }

    int get_vector_varid(const string &name) {
        auto it = vector_vars.find(name);
        if (it == vector_vars.end()) throw runtime_error("Vector variable '" + name + "' not defined");
        return it->second;
    }
};


int main(int argc, char** argv) {
    LibMeshInit init(argc, argv);

    try {
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

        for (size_t i = 0; i < npoints; ++i)
            writer.set_coord(i, points[i]);

        // 3️⃣ Write timestep 0
        writer.set_time(0);
        for (size_t i = 0; i < npoints; ++i) {
            writer.set_scalar("Temperature", i, 300.0+i);
            writer.set_point("Velocity", i, points[i]); // Velocity = point coordinates
            RealTensor T;
            T(0,0)=1.0; T(0,1)=0.1; T(0,2)=0.2;
            T(1,0)=0.1; T(1,1)=1.1; T(1,2)=0.3;
            T(2,0)=0.2; T(2,1)=0.3; T(2,2)=1.2;
            writer.set_tensor("Stress", i, T);
        }

        cout << "NetCDF file 'points_libmesh.nc' written successfully!" << endl;

    } catch (const exception &e) {
        cerr << "Error: " << e.what() << endl;
        return 1;
    }

    return 0;
}

