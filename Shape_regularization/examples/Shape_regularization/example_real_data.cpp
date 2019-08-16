// STL includes.
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <list>
#include <map>

// CGAL includes.
#include <CGAL/property_map.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Shape_regularization.h>
#include <CGAL/linear_least_squares_fitting_2.h>
#include <CGAL/Timer.h>

#include "saver_segments_2.h"

// Typedefs.
typedef CGAL::Exact_predicates_inexact_constructions_kernel Traits;

using Segment_2 = typename Traits::Segment_2;
using Point_2 = typename Traits::Point_2;
using FT = typename Traits::FT;
using Line_2 = typename Traits::Line_2;

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

  void boundary_points_on_line_2(
    const std::vector<Point_2> & item_range,
    const Line_2& line,
    Point_2 &p,
    Point_2 &q) {

    using Vector_2 = typename Traits::Vector_2;
    const FT max_value = FT(1000000000000);

    FT min_proj_value = max_value;
    FT max_proj_value = -max_value;

    const Vector_2 ref_vector = line.to_vector();
    const Point_2& ref_point = item_range[0];
    
    for (std::size_t i = 0; i < item_range.size(); ++i) {
      const Point_2& point = item_range[i];
      
      const Vector_2 curr_vector(ref_point, point);
      const FT value = CGAL::scalar_product(curr_vector, ref_vector);
      
      if (value < min_proj_value) {
        min_proj_value = value;
        p = point; }
      if (value > max_proj_value) {
        max_proj_value = value;
        q = point; }
    }
}

int main(int argc, char *argv[]) {

  std::string path;
  if(argc > 1) {
    path = argv[1];
  }

  CGAL::Timer timer;

  Input_range input_range;

  const std::string data_path = "/media/D/gsoc2019/cgal-dev/Shape_regularization/examples/Shape_regularization/data/real_data_2.xyz";
  std::ifstream file(data_path.c_str(), std::ifstream::in);
  file.precision(15);

  Point_2 p; 
  double stub;
  std::size_t region_index;
  std::map <std::size_t, std::vector<Point_2>> point_map;

  while (!file.eof()) {
    file >> p >> stub >> region_index;
    point_map[region_index].push_back(p);
  }

  std::vector<Line_2> lines;
  for (const auto & m_i : point_map) {
    const std::vector<Point_2> & points = m_i.second;
    Line_2 line;
    linear_least_squares_fitting_2(points.begin(), points.end(), line, CGAL::Dimension_tag<0>());
    lines.push_back(line);
  }

  for (std::size_t i = 0; i < lines.size(); ++i) {
    Point_2 p, q;
    boundary_points_on_line_2(point_map[i], lines[i], p, q) ;
    input_range.push_back(Segment_2(p, q));
  }
  

  std::cout.precision(15);

  std::cout << std::endl;
  std::cout << "BEFORE:" << std::endl;
  for (const auto& segment : input_range)
    std::cout << segment << std::endl;
  std::cout << std::endl;

  Saver saver;
  std::string full_path = path + "example_real_data_before";
  saver.save_segments(input_range, full_path);
 
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

  timer.start();
  shape_regularization_angles.regularize();
  timer.stop();

  std::cout << "Number of modified segments angles: " << regularization_type_angles.number_of_modified_segments()
            << "; Time = " << timer.time() << " sec." << std::endl;
  timer.reset();

  // Regularization for ordinates:
  std::vector <std::vector <std::size_t>> parallel_groups;
  regularization_type_angles.parallel_groups(std::back_inserter(parallel_groups));

  std::cout << "parallel_groups.size() = " << parallel_groups.size() << std::endl;

  const FT bound_ordinates = FT(95) / FT(100);
  Regularization_type_ordinates regularization_type_ordinates(input_range, bound_ordinates);

  neighbor_query.clear();
  for(const auto & group : parallel_groups) {
    neighbor_query.add_group(group);
    regularization_type_ordinates.add_group(group);
  }

  Shape_regularization_ordinates Shape_regularization_ordinates(
    input_range, neighbor_query, regularization_type_ordinates);

  timer.start();
  Shape_regularization_ordinates.regularize();
  timer.stop();
  
  std::cout << "AFTER:" << std::endl;
  for (const auto& segment : input_range)
    std::cout << segment << std::endl;
  std::cout << std::endl;
  
  full_path = path + "example_real_data_after";
  saver.save_segments(input_range, full_path); 

  std::cout << "Number of modified segments ordinates: " << regularization_type_ordinates.number_of_modified_segments()
            << "; Time = " << timer.time() << " sec." << std::endl;


  return EXIT_SUCCESS;
}
