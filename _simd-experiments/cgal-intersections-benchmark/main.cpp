
#include <iostream>

#include <CGAL/point_generators_d.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/point_generators_3.h>

#include <CGAL/Ray_3.h>
#include <CGAL/Triangle_3.h>
#include <CGAL/Bbox_3.h>

using namespace CGAL;
typedef Exact_predicates_inexact_constructions_kernel K;

typedef K::Point_3 Point;
typedef K::Ray_3 Ray;
typedef K::Triangle_3 Triangle;
typedef Bbox_3 Bbox;

typedef Creator_uniform_3<double, Point> Point_creator;
typedef CGAL::Random_points_in_cube_3<Point, Point_creator> Point_generator;

typedef Creator_uniform_2<Point, Ray> Ray_creator;
typedef Creator_uniform_3<Point, Triangle> Triangle_creator;
//typedef Creator_uniform_6<double, Bbox> Bbox_creator;

typedef Join_input_iterator_2<Point_generator, Point_generator, Ray_creator> Ray_generator;
typedef Join_input_iterator_3<Point_generator, Point_generator, Point_generator, Triangle_creator> Triangle_generator;
//typedef Join_input_iterator_2<Point_generator, Point_generator, Bbox_creator> Bbox_generator;


int main() {
  std::size_t Q = 100'000;
  std::size_t T = 100'000;
  std::size_t R = 1'000;

  // All test data will be randomly generated
  // Shapes used for testing will be confined to a cubic region
  CGAL::Random_points_in_cube_3<Point, Point_creator> point_generator(1.0);

  // Generate rays to cast
  std::vector<Ray> ray_queries;
  ray_queries.reserve(Q);
  Ray_generator ray_generator(point_generator, point_generator);
  std::copy_n(ray_generator, Q, std::back_inserter(ray_queries));

  // Generate boxes to cast
  std::vector<Bbox> bbox_queries;
  bbox_queries.reserve(Q);
//  Bbox_generator bbox_generator(point_generator, point_generator);
//  std::copy_n(bbox_generator, Q, std::back_inserter(bbox_queries));

  // Generate boxes to hit
  std::vector<Bbox> bbox_targets;
  bbox_targets.reserve(T);
//  std::copy_n(bbox_generator, T, std::back_inserter(bbox_targets));

  // Generate primitives (triangles) to hit
  std::vector<Triangle> triangle_targets;
  triangle_targets.reserve(T);
  Triangle_generator triangle_generator(point_generator, point_generator, point_generator);
  std::copy_n(triangle_generator, T, std::back_inserter(triangle_targets));

  // Benchmark R times, so that tests are interleaved
  for (int r = 0; r < R; ++r) {

    // Time bbox-bbox intersection

    // Time ray-bbox intersection

    // Time ray-primitive intersection

  }

  // Divide times to produce averages

  // Display results

}