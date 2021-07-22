
#include <iostream>

#include <CGAL/point_generators_d.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/point_generators_3.h>

#include <CGAL/Ray_3.h>
#include <CGAL/Triangle_3.h>
#include <CGAL/Bbox_3.h>

#include <CGAL/Timer.h>

#include "boxed_query.h"

using namespace CGAL;
typedef Exact_predicates_inexact_constructions_kernel K;

typedef K::Point_3 Point;
typedef K::Segment_3 Ray; // FIXME this is a temporary test
typedef K::Triangle_3 Triangle;
typedef Bbox_3 Bbox;

typedef Creator_uniform_3<double, Point> Point_creator;
typedef CGAL::Random_points_in_cube_3<Point, Point_creator> Point_generator;

typedef Creator_uniform_2<Point, Ray> Ray_creator;
typedef Creator_uniform_3<Point, Triangle> Triangle_creator;

struct Bbox_creator {
  Bbox operator()(const Point &min, const Point &size) const {
    return {min.x(), min.y(), min.z(),
            min.x() + size.x(), min.y() + size.y(), min.z() + size.z()};
  }
};

typedef Join_input_iterator_2<Point_generator, Point_generator, Ray_creator> Ray_generator;
typedef Join_input_iterator_3<Point_generator, Point_generator, Point_generator, Triangle_creator> Triangle_generator;
typedef Join_input_iterator_2<Point_generator, Point_generator, Bbox_creator> Bbox_generator;


int main() {
  std::size_t Q = 1'000;
  std::size_t T = 10'000;
  std::size_t R = 100;

  // All test data will be randomly generated
  // Shapes used for testing will be confined to a cubic region

  // Generate rays to cast
  std::vector<Ray> ray_queries;
  ray_queries.reserve(Q);
  Ray_generator ray_generator(Point_generator{}, Point_generator{});
  std::copy_n(ray_generator, Q, std::back_inserter(ray_queries));

  // Generate boxes to cast
  std::vector<Bbox> bbox_queries;
  bbox_queries.reserve(Q);
  Bbox_generator bbox_generator(Point_generator{}, Point_generator{});
  std::copy_n(bbox_generator, Q, std::back_inserter(bbox_queries));

  // Generate boxes to hit
  std::vector<Bbox> bbox_targets;
  bbox_targets.reserve(T);
  std::copy_n(bbox_generator, T, std::back_inserter(bbox_targets));

  // Generate primitives (triangles) to hit
  std::vector<Triangle> primitive_targets;
  primitive_targets.reserve(T);
  Triangle_generator primitive_generator(Point_generator{}, Point_generator{}, Point_generator{});
  std::copy_n(primitive_generator, T, std::back_inserter(primitive_targets));

  // Separate timers for each intersection, keeping track of combined elapsed time of all intersections
  Timer bbox_bbox_timer, ray_bbox_timer, boxed_ray_bbox_timer, ray_primitive_timer;

  // Benchmark R times, so that tests can be interleaved
  for (int r = 0; r < R; ++r) {

    // Draw a progress bar (this benchmark can be pretty slow)
    std::cout << "[";
    for (int i = 0; i < 100; ++i)
      std::cout << (i <= (100 * r / R) ? "=" : " ");
    std::cout << "] (" << r << "/" << R << ")\r" << std::flush;

    // Time bbox-bbox intersection
    bbox_bbox_timer.start();
    for (const auto &query : bbox_queries) {
      for (const auto &target : bbox_targets)
        do_intersect(query, target);
    }
    bbox_bbox_timer.stop();

    // Time ray-bbox intersection
    ray_bbox_timer.start();
    for (const auto &query : ray_queries) {
      for (const auto &target : bbox_targets)
        do_intersect(query, target);
    }
    ray_bbox_timer.stop();

    // Time boxed-ray-bbox intersection
    boxed_ray_bbox_timer.start();
    for (const auto &query : ray_queries) {
      for (const auto &target : bbox_targets)
        do_intersect(Boxed_query<Ray>(query), target);
    }
    boxed_ray_bbox_timer.stop();

    // Time ray-primitive intersection
    ray_primitive_timer.start();
    for (const auto &query : ray_queries) {
      for (const auto &target : primitive_targets)
        do_intersect(query, target);
    }
    ray_primitive_timer.stop();
  }

  // The total count of intersections will be used in calculating averages
  double total_num_intersections = R * Q * T;

  // Display results
  std::cout << "\n"
            << "Intersection times between different types (seconds per intersection)\n"
            << "bbox-bbox: " << bbox_bbox_timer.time() / total_num_intersections << "\n"
            << "ray-bbox: " << ray_bbox_timer.time() / total_num_intersections << "\n"
            << "boxed-ray-bbox: " << boxed_ray_bbox_timer.time() / total_num_intersections << "\n"
            << "ray-primitive: " << ray_primitive_timer.time() / total_num_intersections << "\n"
            << std::endl;

}