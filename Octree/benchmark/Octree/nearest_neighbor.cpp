
#define CGAL_TRACE_STREAM std::cerr

#include <iostream>
#include <fstream>

int main(int argc, char **argv) {

  // Set output file
  std::ofstream file;
  file.open((argc > 1) ? argv[1] : "nearest_neighbor_benchmark.csv");

  file << "Number of Points,Search Time (ms) \n";

  for (size_t num_points = 10; num_points < 10000; num_points *= 1.1) {

    std::cout << num_points << std::endl;

    file << num_points << ",";
    file << 5 << "\n";
  }

  file.close();

  return 0;
}
