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

  const Point_2 a = Point_2(1.0, 3.0);
  const Point_2 b = Point_2(0.918327735082, 2.0221228149508);

  const Point_2 c = Point_2(0.9391957097394, 1.9891280930952);
  const Point_2 d = Point_2(0.9715140165128, 0.983188653235);

  const Point_2 e = Point_2(0.9136318016088, 0.9503266622068);
  const Point_2 f = Point_2(1.996264989341, 1.0264433409682);

  const Point_2 g = Point_2(2.0, 1.0);
  const Point_2 h = Point_2(2.0471841422514, 2.0126669341796);

  const Point_2 i = Point_2(2.0, 2.0);
  const Point_2 j = Point_2(2.0289877554422, 3.0167913534028);

  const Point_2 k = Point_2(0.9865543462772, 2.042315234151);
  const Point_2 l = Point_2(1.0233367445636, 3.0454715510536);

  const Point_2 m = Point_2(2.0369444835817, 2.9916148244169);
  const Point_2 n = Point_2(1.2794349029638, 2.6231845950026);

  const Point_2 o = Point_2(1.637770406629, 2.5930089736413);
  const Point_2 p = Point_2(1.2782107765461, 2.5882986951094);

  const Point_2 q = Point_2(1.2782107765461, 2.331753598516);
  const Point_2 r = Point_2(1.6393664950514, 2.58082650783);

  const Point_2 s = Point_2(1.6069870168406, 2.3417165148885);
  const Point_2 t = Point_2(1.2836564146283, 1.6808881353429);

  const Point_2 u = Point_2(1.638048981006, 1.6933961082739);
  const Point_2 v = Point_2(1.266979117387, 1.3848661093098);
  const Point_2 w = Point_2(1.6463876296266, 1.3806967849995);

  Input_range input_range;
  input_range.push_back(Segment_2(a, b));
  input_range.push_back(Segment_2(c, d));
  input_range.push_back(Segment_2(e, f));
  input_range.push_back(Segment_2(g, h));
  input_range.push_back(Segment_2(i, j));
  input_range.push_back(Segment_2(k, i));
  input_range.push_back(Segment_2(l, m));
  input_range.push_back(Segment_2(n, o));
  input_range.push_back(Segment_2(p, q));
  input_range.push_back(Segment_2(r, s));
  input_range.push_back(Segment_2(q, s));
  input_range.push_back(Segment_2(t, u));
  input_range.push_back(Segment_2(t, v));
  input_range.push_back(Segment_2(v, w));
  input_range.push_back(Segment_2(u, w));

  std::cout.precision(15);

  std::cout << std::endl;
  std::cout << "BEFORE:" << std::endl;
  for (const auto& segment : input_range)
    std::cout << segment << std::endl;
  std::cout << std::endl;

  Saver saver;
  saver.save_segments(input_range, "example_15_segments_before");

  Neighbor_query neighbor_query(input_range);
  std::vector<std::size_t> vec;
  vec.resize(input_range.size());
  std::iota(vec.begin(), vec.end(), 0);
  neighbor_query.add_group(vec);
  
  const FT bound_angles = FT(6.15);
  Regularization_type_angles regularization_type_angles(input_range, bound_angles);
  regularization_type_angles.add_group(vec);

  Shape_regularization_angles shape_regularization_angles(
    input_range, neighbor_query, regularization_type_angles);

  shape_regularization_angles.regularize();

  std::cout << "Number of modified segments angles: " << regularization_type_angles.number_of_modified_segments() << std::endl;

  // Regularization for ordinates:
  std::vector <std::vector <std::size_t>> parallel_groups;
  regularization_type_angles.parallel_groups(std::back_inserter(parallel_groups));

  std::cout << "parallel_groups.size() = " << parallel_groups.size() << std::endl;

  const FT bound_ordinates = FT(0.06);
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
