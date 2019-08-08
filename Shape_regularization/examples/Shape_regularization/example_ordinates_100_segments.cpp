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

#include "saver_segments_2.h"

// Typedefs.
// typedef CGAL::Exact_predicates_inexact_constructions_kernel Traits;
typedef CGAL::Simple_cartesian<double> Traits;
// typedef CGAL::Exact_predicates_exact_constructions_kernel Traits;

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
using Saver = CGAL::Regularization::Saver_segments_2<Traits>;

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

int main() {

  Input_range input_range;
  // input_range.reserve(100);

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

  std::cout.precision(15);

  std::cout << std::endl;
  std::cout << "BEFORE:" << std::endl;
  for (const auto& segment : input_range)
    std::cout << segment << std::endl;
  std::cout << std::endl;

  Saver saver;
  saver.save_segments(input_range, "test_ordinates_100_segments_before");

  const FT tolerance = FT(1);
  Parallel_groups paralle_grouping(input_range, tolerance);

  std::vector <std::vector <std::size_t>> parallel_groups;
  paralle_grouping.parallel_groups(std::back_inserter(parallel_groups));

  std::cout << "parallel_groups.size() = " << parallel_groups.size() << std::endl;

  // Regularization for ordinates:

  const FT bound_ordinates = FT(0.25);
  Regularization_type_ordinates regularization_type_ordinates(input_range, bound_ordinates);

  Neighbor_query neighbor_query(input_range);
  for(const auto & group : parallel_groups) {
    neighbor_query.add_group(group);
    regularization_type_ordinates.add_group(group);
  }

  Shape_regularization_ordinates Shape_regularization_ordinates(
    input_range, neighbor_query, regularization_type_ordinates);
  Shape_regularization_ordinates.regularize();

  std::cout << "Number of modified segments ordinates: " << regularization_type_ordinates.number_of_modified_segments() << std::endl;

  std::cout << "AFTER:" << std::endl;
  for (const auto& segment : input_range)
    std::cout << segment << std::endl;
  std::cout << std::endl;
  saver.save_segments(input_range, "test_ordinates_100_segments_after"); 

  return EXIT_SUCCESS;
}
