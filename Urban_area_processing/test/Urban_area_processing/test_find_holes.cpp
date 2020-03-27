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

#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Triangulation_vertex_base_2.h>

#include <CGAL/Urban_area_processing/internal/property_map.h>
#include <CGAL/Urban_area_processing/utils.h>

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

struct Face_info {

  std::size_t index = std::size_t(-1);
  std::size_t object_index = std::size_t(-1);
  std::size_t label = std::size_t(-1);

  bool used = false;
  bool tagged = false;
};

using Fi = Face_info;
using Fb = CGAL::Triangulation_face_base_with_info_2<Fi, Kernel>;
using Vb = CGAL::Triangulation_vertex_base_2<Kernel>;
    
using Tds = CGAL::Triangulation_data_structure_2<Vb, Fb>;
using Delaunay = CGAL::Delaunay_triangulation_2<Kernel, Tds>;

using Saver = UAP::Saver<Kernel>;

// TODO:
// Fix the API and finish this test.

int main(int argc, char *argv[]) {

  Saver saver;
  const std::string in_path  = argv[1];
  const std::string out_path = "/Users/monet/Documents/gf/urban-area-processing/logs/";
  const FT scale = FT(1) / FT(2); // meters
  const FT noise = FT(1) / FT(4); // meters

  Point_set_3 point_set_3;
  std::vector<std::size_t> building_points;
  UAP::get_building_points(
    in_path, "0", "4", "2", "1", out_path, point_set_3, building_points);
  Point_map point_map(point_set_3, point_set_3.point_map());
  saver.export_points(
    building_points, point_map, out_path + "building_points");

  // finish ...
  Kernel kernel;
  Delaunay delaunay;
  CGAL::Urban_area_processing::find_holes(kernel, delaunay);

  return EXIT_SUCCESS;
}
