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
#include <CGAL/point_generators_2.h>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Triangulation_vertex_base_2.h>

#include <CGAL/Urban_area_processing/utils.h>

#include "include/Saver.h"
#include "include/Utilities.h"

namespace UAP = CGAL::Urban_area_processing;
using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;

using FT = typename Kernel::FT;
using Point_2 = typename Kernel::Point_2;
using Point_3 = typename Kernel::Point_3;

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

using Boundary = std::pair<std::vector<Point_2>, std::size_t>;
using Creator = CGAL::Creator_uniform_2<int, Point_2>;
using Saver = UAP::Saver<Kernel>;

// TODO:
// Improve this test by generating a random triangulation with random holes
// and traversing it.

int main(int argc, char *argv[]) {

  Saver saver;
  const std::string out_path = "/Users/monet/Documents/gf/urban-area-processing/logs/";
  std::cout << std::endl << "* Extracting boundaries..." << std::endl;
  
  std::vector<Point_2> points;
  points.reserve(400);
  CGAL::points_on_square_grid_2(
    255.0, 250, std::back_inserter(points), Creator());

  Delaunay delaunay;
  for (const auto& point : points)
    delaunay.insert(point);

  std::size_t count = 0;
  for (auto face = delaunay.finite_faces_begin();
  face != delaunay.finite_faces_end(); ++face, ++count) {
    face->info().index = count;
    face->info().label = 0;
  }
  saver.export_polygon_soup(delaunay, 0, 
    "/Users/monet/Documents/gf/urban-area-processing/logs/triangulation");

  Kernel kernel;
  std::vector<Boundary> boundaries;
  CGAL::Urban_area_processing::extract_boundary_with_holes_from_triangulation(
    kernel, delaunay, std::back_inserter(boundaries));

  std::cout << "Number of detected boundaries: " << boundaries.size() << std::endl;
  std::cout << std::endl;
  
  for (std::size_t i = 0; i < boundaries.size(); ++i) {
    if (boundaries[i].second == std::size_t(-1))
      saver.export_contour(boundaries[i].first, 
      out_path + "boundary_" + std::to_string(i));
    else
      saver.export_contour(boundaries[i].first, 
      out_path + "hole_" + std::to_string(i));
  }
  return EXIT_SUCCESS;
}
