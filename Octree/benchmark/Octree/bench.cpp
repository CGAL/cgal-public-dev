
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

int main(void) {

  // Create file path from date
  std::stringstream file_path;
  file_path << "../results/";
  auto time_run = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
  file_path << std::put_time(std::localtime(&time_run), "%F"); //ctime(&time_run);
  file_path << ".txt";

  // Create a new file at that path
  std::ofstream results(file_path.str());

  return 0;
}