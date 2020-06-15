
#define CGAL_TRACE_STREAM std::cerr

#include <CGAL/Octree.h>
#include <CGAL/Octree/IO.h>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Point_set_3.h>

#include <fstream>
#include <string>
#include <chrono>
#include <ctime>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_3 Point;
typedef CGAL::Point_set_3<Point> Point_set;
typedef CGAL::Octree
        <Kernel, Point_set, typename Point_set::Point_map, typename Point_set::Vector_map>
        Octree;

#define BUCKET_SIZE 20
#define MAX_DEPTH 10

int main(void) {

  // Create file path from date
  std::stringstream file_path;
  file_path << "../results/";
  auto time_run = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
  file_path << std::put_time(std::localtime(&time_run), "%F"); //ctime(&time_run);
  file_path << ".txt";

  // Create a new file at that path
  std::ofstream file(file_path.str());

  // Add a header to the file
  file << std::put_time(std::localtime(&time_run), "%c") << std::endl;
  file << std::endl;

  // Load the test data
  // TODO

  // Describe the tests
  file << "Benchmark Configuration" << std::endl;
  file << "~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
  file << "Bucket Size: " << BUCKET_SIZE << std::endl;
  file << "Max Depth: " << MAX_DEPTH << std::endl;
  file << "Point Count: " << 5 << std::endl;

  // Show the results
  file << "Results" << std::endl;
  file << "~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
  // TODO

  // Close the file
  file.close();

  return 0;
}