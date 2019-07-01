// STL includes.
#include <vector>
#include <string>
#include <iostream>

// CGAL includes.
#include <CGAL/property_map.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Shape_regularization.h>

// Typedefs.
typedef CGAL::Exact_predicates_inexact_constructions_kernel Traits;

using Segment_2 = typename Traits::Segment_2;
using Point_2 = typename Traits::Point_2;
using FT = typename Traits::FT;

using Input_range = std::vector<Segment_2>;
using Segment_map = CGAL::Identity_property_map<Segment_2>;

using Neighbor_query = CGAL::Regularization::Delaunay_neighbor_query_2<Traits, Input_range, Segment_map>;
using Regularization_type = CGAL::Regularization::Rotated_segments_regularization_2<Traits, Input_range, Segment_map>;
using QP_solver = CGAL::Regularization::Dense_QP_solver<Traits, Input_range>;

using Shape_regularization = CGAL::Regularization::Shape_regularization
  <Traits, Input_range, Neighbor_query, Regularization_type, QP_solver>;

int main() {

  const Point_2 a = Point_2(0.0, 0.0);
  const Point_2 b = Point_2(0.0, 1.0);
  const Point_2 c = Point_2(0.1, 0.0);
  const Point_2 d = Point_2(0.2, 1.0);
  const Point_2 f = Point_2(0.0, 1.1);
  const Point_2 g = Point_2(0.2, 1.1);

  Input_range input_range;
  input_range.push_back(Segment_2(a, b));
  input_range.push_back(Segment_2(c, d));
  input_range.push_back(Segment_2(f, g));

  std::cout << std::endl;
  std::cout << "BEFORE:" << std::endl;
  for (const auto& segment : input_range)
    std::cout << segment << std::endl;
  std::cout << std::endl;

  Segment_map segment_map;

  // Create instances of the classes Neighbor_query and Regularization_type.
  Neighbor_query neighbor_query(input_range, segment_map);
  Regularization_type regularization_type(input_range, segment_map);
  QP_solver qp_solver(input_range);

  Shape_regularization shape_regularization(
    input_range, neighbor_query, regularization_type, qp_solver);
  // Run the algorithm.
  shape_regularization.regularize();
  

  
  std::cout << "AFTER:" << std::endl;
  for (const auto& segment : input_range)
    std::cout << segment << std::endl;
  std::cout << std::endl; 
  
  
  return EXIT_SUCCESS;
}
