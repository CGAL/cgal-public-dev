
#define CGAL_TRACE_STREAM std::cerr

#include <iostream>
#include <CGAL/Octree.h>
#include <CGAL/Octree/IO.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Point_set_3.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/squared_distance_3.h>

#include <chrono>
#include <cassert>

using namespace std::chrono;

typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_3 Point;
typedef Kernel::FT FT;
typedef CGAL::Point_set_3<Point> Point_set;
typedef CGAL::Octree::Octree
        <Point_set, typename Point_set::Point_map>
        Octree;

void naive_vs_accelerated(std::size_t dataset_size) {

  // Create a dataset
  Point_set points;
  CGAL::Random_points_in_cube_3<Point> generator;
  points.reserve(dataset_size);
  for (std::size_t i = 0; i < dataset_size; ++i)
    points.insert(*(generator++));

  // Choose another random point from the same bounds as the dataset
  Point random_point = *(generator++);

  // Use the naive algorithm to find the nearest point in the dataset
  Point naive_nearest = *points.points().begin();
  auto naive_start_time = high_resolution_clock::now();
  {

    FT distance_nearest = std::numeric_limits<FT>::max();
    for (auto &p : points.points()) {

      FT distance_current = CGAL::squared_distance(p, random_point);
      if (distance_current < distance_nearest) {

        distance_nearest = distance_current;
        naive_nearest = p;
      }
    }
  }
  duration<float> naive_elapsed_time = high_resolution_clock::now() - naive_start_time;

  std::cout << "Naive --> "
            << "Closest point to "
            << "(" << random_point << ") "
            << "is "
            << "(" << naive_nearest << ") "
            << "at a distance^2 of "
            << CGAL::squared_distance(naive_nearest, random_point)
            << std::endl;

  // Do the same using the octree
  Point octree_nearest = *generator;
  auto octree_start_time = high_resolution_clock::now();
  {
    // TODO: Write a nearest-neighbor implementation and use it here
  }
  duration<float> octree_elapsed_time = high_resolution_clock::now() - octree_start_time;

  std::cout << "Octree --> "
            << "Closest point to "
            << "(" << random_point << ") "
            << "is "
            << "(" << octree_nearest << ") "
            << "at a distance^2 of "
            << CGAL::squared_distance(octree_nearest, random_point)
            << std::endl;

  // Check that they produce the same answer
  assert(octree_nearest == naive_nearest);

  // Check that the octree was faster
  assert(octree_elapsed_time < naive_elapsed_time);
}

int main(void) {

  naive_vs_accelerated(100);
  naive_vs_accelerated(1000);
  naive_vs_accelerated(10000);
  naive_vs_accelerated(100000);

  return EXIT_SUCCESS;
}
