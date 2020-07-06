
#define CGAL_TRACE_STREAM std::cerr

#include <iostream>
#include <CGAL/Octree.h>
#include <CGAL/Octree/IO.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Point_set_3.h>

#include <cassert>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_3 Point;
typedef CGAL::Point_set_3<Point> Point_set;
typedef CGAL::Octree::Octree
        <Point_set, typename Point_set::Point_map>
        Octree;


void test_1_point() {

  // Define the dataset
  Point_set points;
  points.insert({-1, -1, -1});
  auto point_map = points.point_map();

  // Create the octree
  Octree octree(points, point_map);
  octree.refine(10, 1);

  // Check that the topology matches
  Octree::Node single_node{};
  single_node.value() = octree.root().value();
  assert(single_node == octree.root());

}

void test_2_points() {

  // Define the dataset
  Point_set points;
  points.insert({-1, -1, -1});
  points.insert({1, -1, -1});
  auto point_map = points.point_map();

  // Create the octree
  Octree octree(points, point_map);
  octree.refine(10, 1);

  // The octree should have been split once
  Octree::Node other{};
  other.split();
  assert(other == octree.root());


}

void test_9_points() {

  // TODO
}

int main(void) {


  test_1_point();

  test_2_points();

  test_9_points();

  return EXIT_SUCCESS;
}
