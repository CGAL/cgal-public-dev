// STL includes.
#include <vector>
#include <string>
#include <iostream>

// CGAL includes.
#include <CGAL/property_map.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Shape_regularization/Shape_regularization.h>

// Internal includes.
// #include "include/Saver.h"
// #include "include/Segment_regularizer_2.h"

// Typedefs.
using Traits = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point_2 = typename Traits::Point_2;
using Segment_2 = typename Traits::Segment_2;

// using Saver = CGAL::LOD::Saver<Traits>;
using Identity_map = CGAL::Identity_property_map<Segment_2>;
// using Regularizer = CGAL::LOD::Segment_regularizer_2<Traits>;

int main() {

//   Saver saver;
//   Identity_map identity_map;

  const Point_2 a = Point_2(0.0, 0.0);
  const Point_2 b = Point_2(0.0, 1.0);
  const Point_2 c = Point_2(0.1, 0.0);
  const Point_2 d = Point_2(0.2, 1.0);
  const Point_2 f = Point_2(0.0, 1.1);
  const Point_2 g = Point_2(0.2, 1.1);

  std::vector<Segment_2> segments(3);
  segments[0] = Segment_2(a, b);
  segments[1] = Segment_2(c, d);
  segments[2] = Segment_2(f, g);

  std::cout << std::endl;
  std::cout << "BEFORE:" << std::endl;
  for (const auto& segment : segments)
    std::cout << segment << std::endl;
  std::cout << std::endl;
//   saver.print_segments(segments, identity_map, "before");

//   Regularizer regularizer;
//   regularizer.set_max_angle_in_degrees(10.0);
//   regularizer.set_max_difference_in_meters(0.00001);
//   regularizer.regularize(segments, identity_map);

//   std::cout << std::endl;
  std::cout << "AFTER:" << std::endl;
  for (const auto& segment : segments)
    std::cout << segment << std::endl;
  std::cout << std::endl;
//   saver.print_segments(segments, identity_map, "after");
//   std::cout << std::endl;
  
  return EXIT_SUCCESS;
}
