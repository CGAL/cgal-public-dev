// STL includes.
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <list>

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
using Regularization_type_angles = CGAL::Regularization::Angle_regularization_2<Traits, Input_range, Segment_map>;
using Regularization_type_ordinates = CGAL::Regularization::Ordinate_regularization_2<Traits, Input_range, Segment_map>;

using Shape_regularization_angles = CGAL::Regularization::Shape_regularization
  <Traits, Input_range, Neighbor_query, Regularization_type_angles>;
using Shape_regularization_ordinates = CGAL::Regularization::Shape_regularization
  <Traits, Input_range, Neighbor_query, Regularization_type_ordinates>;
using Parallel_groups = CGAL::Regularization::Parallel_groups_2<Traits, Input_range, Segment_map>;

void print_groups(const std::vector <std::vector <std::size_t>> & groups) {
  std::size_t counter = 0;
  for (const auto & group : groups) {
    std::cout << ++counter << ") ";
    for (const auto & index : group) {
      std::cout << index << " ";
    }
    std::cout << std::endl; 
  }
}

std::stringstream out;

inline std::string data() {
  return out.str();
}

void clear() {
  out.str(std::string());
}

void save(const std::string &file_name, const std::string &extension = ".log") {

  const std::string final_path = file_name + extension;
  std::ofstream file(final_path.c_str(), std::ios_base::out);

  if (!file) std::cerr << std::endl << "ERROR: Error saving log file with the name " << file_name << std::endl << std::endl;

  file << data() << std::endl;
  file.close();
}

template<class Elements, class Segment_map_2>
void save_segments(const Elements &elements, const Segment_map_2 &segment_map_2, const std::string &file_name) {
    
    clear();
    using Const_elements_iterator = typename Elements::const_iterator;

    size_t size = 0;
    for (Const_elements_iterator ce_it = elements.begin(); ce_it != elements.end(); ++ce_it, ++size) {
      const auto &segment = get(segment_map_2, *ce_it);

      out << "v " << segment.source() << " " << 0 << std::endl;
      out << "v " << segment.target() << " " << 0 << std::endl;
      out << "v " << segment.target() << " " << 0 << std::endl;
    }

    for (size_t i = 0; i < size * 3; i += 3)
      out << "f " << i + 1 << " " << i + 2 << " " << i + 3 << std::endl;
    
    save(file_name, ".obj");
}

int main() {

/*
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
// */

  // Test 2.
  // /*
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
  // input_range.erase(input_range.begin() + 16, input_range.begin() + input_range.size());
  // */

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
  /*                                 
  std::cout << "From angle regularization: " << std::endl;
  print_groups(parallel_groups);
  // */
  Parallel_groups parallel_groups_class(input_range);
  std::vector <std::vector <std::size_t>> parallel_groups_2;
  parallel_groups_class.parallel_groups(std::back_inserter(parallel_groups_2));
  /*
  std::cout << "From Parallel_groups: " << std::endl;
  print_groups(parallel_groups_2);
  // */

  Regularization_type_ordinates regularization_type_ordinates(input_range);
  // regularization_type_ordinates.add_groups(parallel_groups);
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

  Segment_map seg_map;
  save_segments(input_range, seg_map, "segment_regularization");
  
  return EXIT_SUCCESS;
}
