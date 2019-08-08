// STL includes.
#include <vector>
#include <string>
#include <iostream>
#include <fstream>

// CGAL includes.
#include <CGAL/property_map.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Shape_regularization.h>

double get_coef_value(const double theta, double & iterator) {
  if (theta == 0 || theta == CGAL_PI / 2 || theta == CGAL_PI || theta == 3 * CGAL_PI / 2) {
    iterator = 0;
  }
  else if (theta == CGAL_PI / 4 || theta == 3 * CGAL_PI / 4 || theta == 5 * CGAL_PI / 4 || theta == 7 * CGAL_PI / 4) {
    iterator = 0.22;
  }
  else if ((theta > 0 && theta < CGAL_PI / 4) || 
           (theta > CGAL_PI / 2 && theta < 3 * CGAL_PI / 4) || 
           (theta > CGAL_PI && theta < 5 * CGAL_PI / 4) || 
           (theta > 3 * CGAL_PI / 2 && theta < 7 * CGAL_PI / 4)) {
             iterator += 0.02;
  }
  else
    iterator -= 0.02;

  if (theta < CGAL_PI) 
    return -1 * iterator;
  return iterator;
}

template<class Traits>
bool test_shape_regularization_segments_2() { 
  using Segment_2 = typename Traits::Segment_2;
  using Point_2 = typename Traits::Point_2;
  using FT = typename Traits::FT;

  using Input_range = std::vector<Segment_2>;
  using Segment_map = CGAL::Identity_property_map<Segment_2>;

  using Neighbor_query = CGAL::Regularization::Delaunay_neighbor_query_2<Traits, Input_range, Segment_map>;
  using Regularization_type_ordinates = CGAL::Regularization::Ordinate_regularization_2<Traits, Input_range, Segment_map>;

  using Shape_regularization_ordinates = CGAL::Regularization::Shape_regularization
      <Traits, Input_range, Neighbor_query, Regularization_type_ordinates>;

  using Parallel_groups = CGAL::Regularization::Parallel_groups_2<Traits, Input_range, Segment_map>;


 Input_range input_range;

  double theta = 0.0;
  double coef = 0.0;
  double iterator = 0.0;
  double theta_step = CGAL_PI / 25.0;

  while(theta < 2 * CGAL_PI) {
    const double st = sin(theta);
    const double ct = cos(theta);
    const Point_2 a = Point_2(0.0, 0.0);
    const Point_2 b = Point_2(ct, st);

    coef = get_coef_value(theta, iterator);
    const Point_2 c = Point_2(ct, st + coef);
    const Point_2 d = Point_2(2 * ct, 2 * st + coef);

    theta += theta_step;

    input_range.push_back(Segment_2(a, b));
    input_range.push_back(Segment_2(c, d));
  }

  assert(input_range.size() == 100);
  if(input_range.size() != 100) {
    return false;
  }

  const FT tolerance = FT(1);
  Parallel_groups paralle_grouping(input_range, tolerance);

  std::vector <std::vector <std::size_t>> parallel_groups;
  paralle_grouping.parallel_groups(std::back_inserter(parallel_groups));

  assert(input_range.size() == 100);
  if(input_range.size() != 100) {;
    return false;
  }
  assert(parallel_groups.size() == 27);
  if(parallel_groups.size() != 27) {
    return false;
  }
  
  const FT bound_ordinates = FT(0.25);
  Regularization_type_ordinates regularization_type_ordinates(input_range, bound_ordinates);

  Neighbor_query neighbor_query(input_range);
  for(const auto & group : parallel_groups) {
    neighbor_query.add_group(group);
    regularization_type_ordinates.add_group(group);
  }

  Shape_regularization_ordinates shape_regularization_ordinates(
    input_range, neighbor_query, regularization_type_ordinates);
  shape_regularization_ordinates.regularize();

  assert(input_range.size() == 100);
  if(input_range.size() != 100) {;
    return false;
  }
  
  const std::size_t modified_seg_ord = regularization_type_ordinates.number_of_modified_segments();
  assert(modified_seg_ord == 100);
  if(modified_seg_ord != 100) {
    return false;
  }

  return true;
}

int main() {
  bool cartesian_double_test_success = true;
  if (!test_shape_regularization_segments_2< CGAL::Simple_cartesian<double> >()) 
    cartesian_double_test_success = false;
      
  std::cout << "cartesian_double_test_success: " << cartesian_double_test_success << std::endl;
  assert(cartesian_double_test_success);

  const bool success = cartesian_double_test_success;
  return (success) ? EXIT_SUCCESS : EXIT_FAILURE;
}
