// STL includes.
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <cmath>

// CGAL includes.
#include <CGAL/property_map.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Shape_regularization.h>


template<class Traits>
bool test_shape_regularization_segments_2() { 
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

  Input_range input_range;

  const std::string testpath = "/media/D/gsoc2019/cgal-dev/Shape_regularization/test/Shape_regularization/data/segments_216.polylines";
  // std::cout << testpath << std::endl;
  std::ifstream file(testpath.c_str(), std::ifstream::in);
  file.precision(15);

  Point_2 s, t; double stub;
  while (!file.eof()) {
    file >> stub >> s >> stub >> t >> stub;
    input_range.push_back(Segment_2(s, t));
  }
  input_range.erase(input_range.begin() + input_range.size() - 1);


  assert(input_range.size() == 216);
  if(input_range.size() != 216) {
    return false;
  }

  Neighbor_query neighbor_query(input_range);
  std::vector<std::size_t> vec;
  vec.resize(input_range.size());
  std::iota(vec.begin(), vec.end(), 0);
  neighbor_query.add_group(vec);
  
  const FT bound_angles = FT(60);
  Regularization_type_angles regularization_type_angles(input_range, bound_angles);
  regularization_type_angles.add_group(vec);

  Shape_regularization_angles shape_regularization_angles(
    input_range, neighbor_query, regularization_type_angles);

  shape_regularization_angles.regularize();

  std::vector <std::vector <std::size_t>> parallel_groups;
  regularization_type_angles.parallel_groups(std::back_inserter(parallel_groups));

  assert(input_range.size() == 216);
  if(input_range.size() != 216) {;
    return false;
  }
  assert(parallel_groups.size() == 30);
  if(parallel_groups.size() != 30) {
    return false;
  }
  
  const std::size_t modified_seg_ang = regularization_type_angles.number_of_modified_segments();
  assert(modified_seg_ang == 216);
  if(modified_seg_ang != 216) {
    return false;
  }

  // Regularization for ordinates:

  const FT bound_ordinates = FT(5);
  Regularization_type_ordinates regularization_type_ordinates(input_range, bound_ordinates);

  neighbor_query.clear();
  for(const auto & group : parallel_groups) {
    neighbor_query.add_group(group);
    regularization_type_ordinates.add_group(group);
  }

  Shape_regularization_ordinates shape_regularization_ordinates(
    input_range, neighbor_query, regularization_type_ordinates);
  shape_regularization_ordinates.regularize();

  assert(input_range.size() == 216);
  if(input_range.size() != 216) {
    return false;
  }

  const std::size_t modified_seg_ord = regularization_type_ordinates.number_of_modified_segments();
  assert(modified_seg_ord == 202);
  if(modified_seg_ord != 202) {
    return false;
  }

  return true;
}

int main() {
  bool exact_exact_test_success = true;
  if (!test_shape_regularization_segments_2<CGAL::Exact_predicates_exact_constructions_kernel>()) 
    exact_exact_test_success = false;
    
  std::cout << "exact_exact_test_success: " << exact_exact_test_success << std::endl;
  assert(exact_exact_test_success);

  const bool success = exact_exact_test_success;
  return (success) ? EXIT_SUCCESS : EXIT_FAILURE;
}
