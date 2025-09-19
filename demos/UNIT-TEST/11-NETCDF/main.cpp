
#include "util/NetCDFWriter.h"

using namespace util;

/**
 *
 *
 *    MAIN
 *
 */

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

    writer.define_variable(TEMPERATURE);
    writer.define_variable(PRESSURE);
    writer.define_variable(VELOCITY);
    writer.define_variable(STRESS);

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
      writer.set_curr_time(time);

      // Each processor writes data for its own local points
      for (uint i = start_idx; i < end_idx; ++i) 
      {
        writer.set_curr_pt(i);
        double temp = generate_temperature(i, time);
        writer.set_value(TEMPERATURE, temp);

        double pressure = 101325.0 + temp * 100; // Pressure depends on temperature
        writer.set_value(PRESSURE, pressure);

        // Set vector values
        Point velocity = generate_velocity(i, time);
        writer.set_value(VELOCITY, velocity);

        // Set tensor values
        RealTensorValue stress = generate_stress(i, time);
        writer.set_value(STRESS, stress);
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
