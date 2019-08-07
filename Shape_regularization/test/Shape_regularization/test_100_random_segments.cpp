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

// for random generator
#include <algorithm>
#include <CGAL/point_generators_2.h>
#include <CGAL/function_objects.h>
#include <CGAL/Join_input_iterator.h>
#include <CGAL/Counting_iterator.h>


template<class Traits>
bool test_shape_regularization_segments_2() { 
  using Segment_2 = typename Traits::Segment_2;
  using Point_2 = typename Traits::Point_2;
  using FT = typename Traits::FT;

  using Input_range = std::vector<Segment_2>;
  using Segment_map = CGAL::Identity_property_map<Segment_2>;

  using Neighbor_query = CGAL::Regularization::Delaunay_neighbor_query_2<Traits, Input_range, Segment_map>;
  using Regularization_type_angles = CGAL::Regularization::Angle_regularization_2<Traits, Input_range, Segment_map>;
  using Shape_regularization_angles = CGAL::Regularization::Shape_regularization
    <Traits, Input_range, Neighbor_query, Regularization_type_angles>;

  typedef CGAL::Points_on_segment_2<Point_2>              PG;
  typedef CGAL::Creator_uniform_2< Point_2, Segment_2>    Creator;
  typedef CGAL::Join_input_iterator_2< PG, PG, Creator>   Segm_iterator;
  typedef CGAL::Counting_iterator<Segm_iterator, Segment_2>  Count_iterator;

  Input_range input_range;
  input_range.reserve(100);

  // A horizontal like fan.
  PG p1( Point_2(-250, -50), Point_2(-250, 50),50);   // Point generator.
  PG p2( Point_2( 250,-250), Point_2( 250,250),50);
  Segm_iterator  t1( p1, p2);                     // Segment generator.
  Count_iterator t1_begin( t1);                   // Finite range.
  Count_iterator t1_end( t1, 50);
  std::copy( t1_begin, t1_end, std::back_inserter(input_range));
  // A vertical like fan.
  PG p3( Point_2( -50,-250), Point_2(  50,-250),50);
  PG p4( Point_2(-250, 250), Point_2( 250, 250),50);
  Segm_iterator  t2( p3, p4);
  Count_iterator t2_begin( t2);
  Count_iterator t2_end( t2, 50);
  std::copy( t2_begin, t2_end, std::back_inserter(input_range));

  assert(input_range.size() == 100);
  if(input_range.size() != 100) {
    return false;
  }

  Neighbor_query neighbor_query(input_range);
  std::vector<std::size_t> vec;
  vec.resize(input_range.size());
  std::iota(vec.begin(), vec.end(), 0);
  neighbor_query.add_group(vec);
  
  const FT bound_angles = FT(39);
  Regularization_type_angles regularization_type_angles(input_range, bound_angles);
  regularization_type_angles.add_group(vec);

  Shape_regularization_angles shape_regularization_angles(
    input_range, neighbor_query, regularization_type_angles);

  shape_regularization_angles.regularize();

  std::vector <std::vector <std::size_t>> parallel_groups;
  regularization_type_angles.parallel_groups(std::back_inserter(parallel_groups));

  assert(input_range.size() == 100);
  if(input_range.size() != 100) {;
    return false;
  }
  assert(parallel_groups.size() == 2);
  if(parallel_groups.size() != 2) {
    return false;
  }
  
  const std::size_t modified_seg_ang = regularization_type_angles.number_of_modified_segments();
  assert(modified_seg_ang == 100);
  if(modified_seg_ang != 100) {
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

  // ------>

  bool exact_inexact_test_success = true;
  if (!test_shape_regularization_segments_2<CGAL::Exact_predicates_inexact_constructions_kernel>()) 
    exact_inexact_test_success = false;
    
  std::cout << "exact_inexact_test_success: " << exact_inexact_test_success << std::endl;
  assert(exact_inexact_test_success);

  // ------>

  bool exact_exact_test_success = true;
  if (!test_shape_regularization_segments_2<CGAL::Exact_predicates_exact_constructions_kernel>()) 
    exact_exact_test_success = false;
    
  std::cout << "exact_exact_test_success: " << exact_exact_test_success << std::endl;
  assert(exact_exact_test_success);

  const bool success = cartesian_double_test_success && exact_inexact_test_success && exact_exact_test_success;
  return (success) ? EXIT_SUCCESS : EXIT_FAILURE;
}
