#include <netcdf>
#include <iostream>

using namespace netCDF;

int main() {
    try {
        // Create a classic NetCDF file
        NcFile file("simple_test.nc", NcFile::replace, NcFile::classic);

        // Define a single-dimension variable (length 1)
        NcDim dim = file.addDim("scalar_dim", 1);
        NcVar var = file.addVar("my_value", ncFloat, dim);

        // Write a single value
        float value = 42.0f;
        var.putVar(&value);

        // Read it back
        float readValue;
        file.getVar("my_value").getVar(&readValue);

        std::cout << "Value read from NetCDF file: " << readValue << "\n";

    } catch (netCDF::exceptions::NcException &e) {
        std::cerr << "NetCDF error: " << e.what() << "\n";
        return 1;
    }

    return 0;
}

