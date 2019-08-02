// STL includes.
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <list>

// CGAL includes.
#include <CGAL/property_map.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Shape_regularization.h>

// Typedefs.
typedef CGAL::Exact_predicates_inexact_constructions_kernel Traits;
// typedef CGAL::Simple_cartesian<double> Traits;
// typedef CGAL::Exact_predicates_exact_constructions_kernel Traits;

using Segment_2 = typename Traits::Segment_2;
using Point_2 = typename Traits::Point_2;
using FT = typename Traits::FT;

using Input_range = std::vector<Segment_2>;
using Segment_map = CGAL::Identity_property_map<Segment_2>;

using Neighbor_query = CGAL::Regularization::Delaunay_neighbor_query_2<Traits, Input_range, Segment_map>;
using Regularization_type_angles = CGAL::Regularization::Angle_regularization_2<Traits, Input_range, Segment_map>;
using Regularization_type_ordinates = CGAL::Regularization::Ordinate_regularization_2<Traits, Input_range, Segment_map>;

using Shape_regularization_angles = CGAL::Regularization::Shape_regularization
  <Traits, Input_range, Neighbor_query, Regularization_type_angles>;
using Shape_regularization_ordinates = CGAL::Regularization::Shape_regularization
  <Traits, Input_range, Neighbor_query, Regularization_type_ordinates>;
using Parallel_groups = CGAL::Regularization::Parallel_groups_2<Traits, Input_range, Segment_map>;

int main() {

  Input_range input_range;
  const std::string testpath = "/media/D/gsoc2019/cgal-dev/Shape_regularization/examples/Shape_regularization/data/test.polylines";
  std::cout << testpath << std::endl;
  std::ifstream file(testpath.c_str(), std::ifstream::in);
  file.precision(15);

  Point_2 s, t; double stub;
  while (!file.eof()) {
    file >> stub >> s >> stub >> t >> stub;
    input_range.push_back(Segment_2(s, t));
  }
  input_range.erase(input_range.begin() + input_range.size() - 1);

  std::cout.precision(15);

  std::cout << std::endl;
  std::cout << "BEFORE:" << std::endl;
  for (const auto& segment : input_range)
    std::cout << segment << std::endl;
  std::cout << std::endl;

  // Create instances of the classes Neighbor_query and Regularization_type.
  // Neighbor_query neighbor_query_angles(input_range);
  Neighbor_query neighbor_query(input_range);
  Regularization_type_angles regularization_type_angles(input_range);

  Shape_regularization_angles shape_regularization_angles(
    input_range, neighbor_query, regularization_type_angles);

  shape_regularization_angles.regularize();


  // Regularization for ordinates:
  std::vector <std::vector <std::size_t>> parallel_groups;
  regularization_type_angles.parallel_groups(std::back_inserter(parallel_groups));

  Regularization_type_ordinates regularization_type_ordinates(input_range);

  neighbor_query.clear();
  for(const auto & group : parallel_groups) {
    neighbor_query.add_group(group);
    regularization_type_ordinates.add_group(group);
  }

  Shape_regularization_ordinates Shape_regularization_ordinates(
    input_range, neighbor_query, regularization_type_ordinates);

  Shape_regularization_ordinates.regularize();
  
  std::cout << "AFTER:" << std::endl;
  for (const auto& segment : input_range)
    std::cout << segment << std::endl;
  std::cout << std::endl; 
  
  return EXIT_SUCCESS;
}
