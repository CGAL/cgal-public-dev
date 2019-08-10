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

#include "saver_segments_2.h"

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

using Saver = CGAL::Regularization::Saver_segments_2<Traits>;


int main() {

  const Point_2 a = Point_2(1.0, 1.0);
  const Point_2 a1 = Point_2(1.968759150567688, 1.599174332100224);
  const Point_2 a2 = Point_2(0.93, 1.0);
  const Point_2 b = Point_2(0.925377816338188, 2.995179914344531);
  const Point_2 c = Point_2(1.0, 3.0);
  const Point_2 d = Point_2(1.066662126401646, 4.951894853937938);
  const Point_2 e = Point_2(1.0, 5.0);
  const Point_2 f = Point_2(2.95, 4.930389132419256);
  const Point_2 g = Point_2(3.0, 4.95);
  const Point_2 h = Point_2(2.934996832312081, 3.008203183890471);
  const Point_2 i = Point_2(3.085452511148867, 3.003266787827703);
  const Point_2 j = Point_2(2.96978240203571, 1.002004749629305);
  const Point_2 j1 = Point_2(2.86, 1.002004749629305);
  const Point_2 k = Point_2(0.948866110654676, 3.033161728487555);
  const Point_2 l = Point_2(2.9, 3.0);
  const Point_2 m = Point_2(1.6, 4.0);
  const Point_2 n = Point_2(1.932136680786834, 4.36471855871896);
  const Point_2 o = Point_2(1.598613201104051, 3.982686390724744);
  const Point_2 p = Point_2(2.018220171854482, 3.686595362285412);
  const Point_2 q = Point_2(1.951872279478803, 4.363094734768067);
  const Point_2 r = Point_2(2.290848420455327, 4.0541544543844);
  const Point_2 s = Point_2(2.30451735553157, 4.045054694344393);
  const Point_2 t = Point_2(1.642059717842882, 1.928505579230186);
  const Point_2 u = Point_2(1.993860111389907, 2.247986994205749);
  const Point_2 w = Point_2(2.259099631673991, 1.919966912585693);
  const Point_2 z = Point_2(1.62984590502084, 1.923077217975945);

  std::vector<std::size_t> group1; 
  std::vector<std::size_t> group2; // the top romb
  std::vector<std::size_t> group3; // the bottom romb

  Input_range input_range;
  input_range.push_back(Segment_2(a, b));
  input_range.push_back(Segment_2(c, d));
  input_range.push_back(Segment_2(e, f));
  input_range.push_back(Segment_2(g, h));
  input_range.push_back(Segment_2(i, j));
  input_range.push_back(Segment_2(k, l));
  input_range.push_back(Segment_2(a2, j1));
  group1.push_back(0);
  group1.push_back(1);
  group1.push_back(2);
  group1.push_back(3);
  group1.push_back(4);
  group1.push_back(5);
  group1.push_back(6);

  input_range.push_back(Segment_2(m, n));
  input_range.push_back(Segment_2(o, p));
  input_range.push_back(Segment_2(q, r));
  input_range.push_back(Segment_2(p, s));
  group2.push_back(7);
  group2.push_back(8);
  group2.push_back(9);
  group2.push_back(10);

  input_range.push_back(Segment_2(t, u));
  input_range.push_back(Segment_2(u, w));
  input_range.push_back(Segment_2(z, a1));
  input_range.push_back(Segment_2(w, a1));
  group3.push_back(11);
  group3.push_back(12);
  group3.push_back(13);
  group3.push_back(14);

  std::cout.precision(15);

  std::cout << std::endl;
  std::cout << "BEFORE:" << std::endl;
  for (const auto& segment : input_range)
    std::cout << segment << std::endl;
  std::cout << std::endl;

  Saver saver;
  saver.save_segments(input_range, "example_15_segments_before");

  Neighbor_query neighbor_query(input_range);
  neighbor_query.add_group(group1);
  neighbor_query.add_group(group2);
  neighbor_query.add_group(group3);
  
  const FT bound_angles = FT(385) / FT(100);
  Regularization_type_angles regularization_type_angles(input_range, bound_angles);
  regularization_type_angles.add_group(group1);
  regularization_type_angles.add_group(group2);
  regularization_type_angles.add_group(group3);

  Shape_regularization_angles shape_regularization_angles(
    input_range, neighbor_query, regularization_type_angles);

  shape_regularization_angles.regularize();

  std::cout << "Number of modified segments angles: " << regularization_type_angles.number_of_modified_segments() << std::endl;

  // Regularization for ordinates:
  std::vector <std::vector <std::size_t>> parallel_groups;
  regularization_type_angles.parallel_groups(std::back_inserter(parallel_groups));

  std::cout << "parallel_groups.size() = " << parallel_groups.size() << std::endl;

  const FT bound_ordinates = FT(1) / FT(10);
  Regularization_type_ordinates regularization_type_ordinates(input_range, bound_ordinates);

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
  saver.save_segments(input_range, "example_15_segments_after"); 

  std::cout << "Number of modified segments ordinates: " << regularization_type_ordinates.number_of_modified_segments() << std::endl;


  return EXIT_SUCCESS;
}
