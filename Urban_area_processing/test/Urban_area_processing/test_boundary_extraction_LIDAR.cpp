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

#include <CGAL/Urban_area_processing/property_map.h>
#include <CGAL/Urban_area_processing/Boundary_extraction.h>

#include "include/Saver.h"
#include "include/Utilities.h"

namespace UAP = CGAL::Urban_area_processing;
using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;

using FT = typename Kernel::FT;
using Point_2 = typename Kernel::Point_2;
using Point_3 = typename Kernel::Point_3;

using Point_set_3 = CGAL::Point_set_3<Point_3>;
using Point_map_3 = typename Point_set_3::Point_map;

using Point_map = UAP::Item_property_map<Point_set_3, Point_map_3>;
using Boundary_extraction_LIDAR = UAP::Boundary_extraction_LIDAR<
  Kernel, std::vector<std::size_t>, Point_map>;

using Saver = UAP::Saver<Kernel>;

int main(int argc, char *argv[]) {

  Saver saver;
  const std::string in_path  = argv[1];
  const std::string out_path = "/Users/monet/Documents/gf/urban-area-processing/logs/";

  Point_set_3 point_set_3;
  std::vector<std::size_t> building_points;
  UAP::get_building_points(
    in_path, "0", "1", "2", "3", out_path, point_set_3, building_points);
  Point_map point_map(point_set_3, point_set_3.point_map());
  saver.export_points(
    building_points, point_map, out_path + "building_points");

  Boundary_extraction_LIDAR extractor(
    building_points, point_map);

  std::vector<Point_2> boundary;
  extractor.extract(std::back_inserter(boundary));
  saver.export_contour(boundary, out_path + "boundary_LIDAR");

  return EXIT_SUCCESS;
}
