
#define CGAL_TRACE_STREAM std::cerr

#include <iostream>
#include <fstream>
#include <chrono>

#include "util.h"

using std::chrono::high_resolution_clock;
using std::chrono::duration_cast;
using std::chrono::microseconds;

int main(int argc, char **argv) {

  // Set output file
  std::ofstream file;
  file.open((argc > 1) ? argv[1] : "construction_benchmark.csv");

  // Add header for CSV
  file << "Number of Points,Build Time (ms) \n";

  // Perform tests for various dataset sizes
  for (size_t num_points = 10; num_points < 10000; num_points *= 1.1) {

    auto start = high_resolution_clock::now();

    // TODO

    auto end = high_resolution_clock::now();

    file << num_points << ",";
    file << duration_cast<microseconds>(end - start).count() << "\n";
  }

  file.close();

  return 0;
}