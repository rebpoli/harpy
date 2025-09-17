// save as write_netcdf.cpp
#include <vector>
#include <string>
#include <iostream>

#include <netcdf.h>

using namespace netCDF;
using namespace netCDF::exceptions;

int main() {
    try {
        // Create a new NetCDF4 file (overwrite if exists)
        NcFile dataFile("data.nc", NcFile::replace, NcFile::nc4);

        // Define dimensions
        int N = 10;
        NcDim dim = dataFile.addDim("x", N);

        // Define variables
        NcVar var_values = dataFile.addVar("values", ncDouble, dim);
        NcVar var_names = dataFile.addVar("names", ncChar, dim);

        // Write values
        std::vector<double> values(N);
        std::vector<std::string> names = {"a","b","c","d","e","f","g","h","i","j"};
        for (int i = 0; i < N; i++) values[i] = i * 1.5;

        var_values.putVar(values.data());

        // Write names as fixed-length strings
        for (int i = 0; i < N; i++) {
            var_names.putVar({i}, {1}, names[i].c_str());
        }

        std::cout << "NetCDF4 file written successfully!\n";
    }
    catch (NcException& e) {
        e.what();
        return 1;
    }
    return 0;
}

