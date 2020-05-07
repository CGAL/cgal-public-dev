// STL includes.
#include <string>
#include <vector>
#include <utility>
#include <fstream>
#include <iostream>
#include <iterator>
#include <cassert>

// CGAL includes.
#include <CGAL/assertions.h>
#include <CGAL/property_map.h>
#include <CGAL/Point_set_3.h>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <CGAL/Urban_area_processing/internal/property_map.h>
#include <CGAL/Urban_area_processing/Boundaries.h>

#include "include/Saver.h"
#include "include/Utilities.h"

namespace UAP = CGAL::Urban_area_processing;
using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;

using FT = typename Kernel::FT;
using Point_2 = typename Kernel::Point_2;
using Point_3 = typename Kernel::Point_3;

using Point_set_3 = CGAL::Point_set_3<Point_3>;
using Point_map_3 = typename Point_set_3::Point_map;

using Point_map = UAP::internal::Item_property_map<Point_set_3, Point_map_3>;
using Boundary_extraction = UAP::Boundary_extraction_with_repair<
  Kernel, std::vector<std::size_t>, Point_map>;

using Boundary = std::pair<std::vector<Point_2>, std::size_t>;
using Saver = UAP::Saver<Kernel>;

// TODO:
// I think, I should merge this test later with the LIDAR version test.

int main(int argc, char *argv[]) {

  Saver saver;
  const std::string in_path  = argv[1];
  const std::string out_path = "/Users/monet/Documents/gf/urban-area-processing/logs/";
  const FT scale = FT(1); // meters
  const FT noise = FT(1); // meters
  const FT min_length_2 = FT(2);  // meters
  const FT max_angle_2  = FT(25); // degrees
  const FT max_angle_3  = FT(25); // degrees

  Point_set_3 point_set_3;
  std::vector<std::size_t> building_points;
  UAP::get_building_points(
    in_path, "0", "4", "2", "1", out_path, point_set_3, building_points);
  Point_map point_map(point_set_3, point_set_3.point_map());
  saver.export_points(
    building_points, point_map, out_path + "building_points");

  Boundary_extraction extractor(
    building_points, point_map, 
    scale, noise, min_length_2, max_angle_2, max_angle_3);

  std::vector<Boundary> boundaries;
  extractor.extract(std::back_inserter(boundaries));

  std::cout << "Number of detected boundaries: " << boundaries.size() << std::endl;
  std::cout << std::endl;
  
  for (std::size_t i = 0; i < boundaries.size(); ++i) {
    if (boundaries[i].second == std::size_t(-1))
      saver.export_contour(boundaries[i].first, 
      out_path + "boundary_with_repair_" + std::to_string(i));
    else
      saver.export_contour(boundaries[i].first, 
      out_path + "hole_with_repair_" + std::to_string(i));
  }
  return EXIT_SUCCESS;
}
