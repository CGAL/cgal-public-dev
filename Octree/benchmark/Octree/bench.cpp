
#define CGAL_TRACE_STREAM std::cerr

#include <CGAL/Octree.h>
#include <CGAL/Octree/IO.h>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Point_set_3.h>
#include <CGAL/Point_set_3/IO.h>

#include <fstream>
#include <string>
#include <chrono>
#include <ctime>
#include <string>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_3 Point;
typedef CGAL::Point_set_3<Point> Point_set;
typedef CGAL::Octree
        <Kernel, Point_set, typename Point_set::Point_map>
        Octree;

#define DATA_FILE_PATH "../data/archer_statue_scan.ply"
#define BUCKET_SIZE 20
#define MAX_DEPTH 10
#define NUM_RUNS 40

std::chrono::duration<double> benchmark(Point_set points) {

  std::cout << "Benchmarking" << std::endl;

  std::chrono::duration<double> elapsed(0);

  for (int i = 0; i < NUM_RUNS; ++i) {

    auto points_copy = points;
    auto point_map = points_copy.point_map();

    auto start = std::chrono::high_resolution_clock::now();

    Octree octree(points, point_map);
    octree.refine(MAX_DEPTH, BUCKET_SIZE);

    auto end = std::chrono::high_resolution_clock::now();

    std::cout << "." << std::flush;

    elapsed += (end - start) / NUM_RUNS;
  }


  std::cout << std::endl << "Averaged "
            << std::chrono::duration_cast<std::chrono::milliseconds>(elapsed).count() << " ms" << std::endl;
  return elapsed;
}

int main(void) {

  // Create file path from date
  std::stringstream file_path;
  file_path << "../results/";
  auto time_run = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
  file_path << std::put_time(std::localtime(&time_run), "%F"); //ctime(&time_run);
  file_path << ".txt";

  // Create a new file at that path
  std::ofstream file(file_path.str());

  // Add a header to the file with the time it was run
  file << std::put_time(std::localtime(&time_run), "%c") << std::endl;
  file << std::endl;

  // Load the test data
  std::ifstream stream(DATA_FILE_PATH);
  Point_set points;
  stream >> points;

  // Describe the tests
  file << "Benchmark Configuration" << std::endl;
  file << "~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
  file << "Bucket Size: " << BUCKET_SIZE << std::endl;
  file << "Max Depth: " << MAX_DEPTH << std::endl;
  file << "Dataset Source: " << DATA_FILE_PATH << std::endl;
  file << "Point Count: " << points.number_of_points() << std::endl;
  file << "Number of Runs: " << NUM_RUNS << std::endl;
  file << std::endl;

  // Show the results
  file << "Results" << std::endl;
  file << "~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
  file << "Tree construction time: "
       << std::chrono::duration_cast<std::chrono::milliseconds>(benchmark(points)).count() << " ms" << std::endl;
  file << std::endl;

  // Leave room for notes
  file << "Notes" << std::endl;
  file << "~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
  file << std::endl;

  // Close the file
  file.close();

  return 0;
}