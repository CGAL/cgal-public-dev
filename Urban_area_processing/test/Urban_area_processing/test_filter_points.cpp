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

#include <CGAL/Urban_area_processing/utils.h>
#include <CGAL/Urban_area_processing/internal/property_map.h>
#include <CGAL/Urban_area_processing/internal/Shape_detection/Estimate_normals_3.h>
#include <CGAL/Urban_area_processing/internal/Shape_detection/Sphere_neighbor_query.h>
#include <CGAL/Urban_area_processing/internal/Tools/Extract_vertical_points_3.h>

#include "include/Saver.h"
#include "include/Utilities.h"

namespace UAP = CGAL::Urban_area_processing;
using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;

using FT = typename Kernel::FT;
using Point_2 = typename Kernel::Point_2;
using Point_3 = typename Kernel::Point_3;
using Vector_3 = typename Kernel::Vector_3;

using Point_set_3 = CGAL::Point_set_3<Point_3>;
using Point_map_3 = typename Point_set_3::Point_map;

using Indices = std::vector<std::size_t>;
using Point_map = UAP::internal::Item_property_map<Point_set_3, Point_map_3>;

using Sphere_neighbor_query = UAP::internal::Sphere_neighbor_query<
  Kernel, Indices, Point_map>;
using Estimate_normals_3 = UAP::internal::Estimate_normals_3<
  Kernel, Indices, Sphere_neighbor_query>;
using Vertical_condition = 
  UAP::internal::Extract_vertical_points_3<Kernel>;

using Saver = UAP::Saver<Kernel>;
using Color = CGAL::Color;

void estimate_normals(
  const Indices& building_points,
  const Point_map point_map,
  const FT scale,
  std::vector<Vector_3>& normals) {

  Sphere_neighbor_query neighbor_query(
    building_points, scale, point_map);
  Estimate_normals_3 estimator(
    building_points, neighbor_query);
  estimator.get_normals(normals);
  assert(normals.size() == building_points.size());
}

// TODO:
// Simplify this example by using already implemented functions from other
// CGAL packages instead of my internal functions/classes.

int main(int argc, char *argv[]) {

  Saver saver;
  const std::string in_path  = argv[1];
  const std::string out_path = "/Users/monet/Documents/gf/urban-area-processing/logs/";
  const FT scale = FT(1) / FT(2); // meters
  const FT max_angle_3 = FT(25); // degrees

  Point_set_3 point_set_3;
  std::vector<std::size_t> building_points;
  UAP::get_building_points(
    in_path, "0", "4", "2", "1", out_path, point_set_3, building_points);
  Point_map point_map(point_set_3, point_set_3.point_map());
  saver.export_points(
    building_points, point_map, out_path + "building_points");

  std::vector<Vector_3> normals;
  estimate_normals(
    building_points, point_map, scale, normals);
  const Vertical_condition vertical_condition(
    normals, max_angle_3);

  Kernel kernel;
  std::vector<Point_3> points;
  UAP::filter_points(
    kernel, building_points, vertical_condition, point_map, 
    std::back_inserter(points));

  std::cout << "Number of extracted points: " << points.size() << std::endl;
  std::cout << std::endl;
  
  saver.export_points(points, Color(0, 0, 0), out_path + "vertical_points");
  return EXIT_SUCCESS;
}
