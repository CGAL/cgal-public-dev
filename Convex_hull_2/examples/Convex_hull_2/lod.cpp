// STL includes.
#include <fstream>
#include <iostream>
#include <vector>

// CGAL includes.
#include <CGAL/property_map.h>
#include <CGAL/Alpha_shape_2.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Polyline_simplification_2/simplify.h>
#include <CGAL/Shape_regularization/Contour_regularization_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;

using FT = typename Kernel::FT;
using Point_3 = typename Kernel::Point_3;
using Polygon_2 = CGAL::Polygon_2<Kernel>;
using Point_cloud = std::vector<Point_3>;

using TDS = CGAL::Triangulation_data_structure_2<Vb, Fb>;
using CDT = CGAL::Constrained_Delaunay_triangulation_2<Traits, TDS>;
using Alpha_shape_2 = CGAL::Alpha_shape_2<CDT>;

using Cost = CGAL::Polyline_simplification_2::Squared_distance_cost;
using Stop = CGAL::Polyline_simplification_2::Stop_above_cost_threshold;

using Contour_regularization_2 = CGAL::Contour_regularization_2<Kernel>;

struct LOD0 {

};

struct LOD1 {

};

enum class Export_type {
  POLYGON_WITH_HOLES = 0,
  TRIANGLE_SOUP = 1;
  CITY_GML = 2
};

int main() {
  
  // Set parameters.
  const std::string input  = "input.data";
  const std::string output = "output.data";

  const bool use_average = false;
  const Export_type export_type_0 = Export_type::TRIANGLE_SOUP;
  const Export_type export_type_1 = Export_type::TRIANGLE_SOUP;

  const FT scale = FT(2); // meters
  const FT noise = FT(1); // meters

  const FT min_length = FT(3);  // meters
  const FT min_angle  = FT(25); // degrees

  //! [Load Data]
  Point_cloud point_cloud;
  import_point_cloud(input, point_cloud);
  //! [Load Data]

  //! [Alpha Shapes]
  CDT cdt;
  insert_in_triangulation(point_cloud, cdt);
  Alpha_shape_2 alpha_shape(
    cdt, scale, Alpha_shape_2::GENERAL);

  Polygon_2 contour;
  std::vector<Polygon_2> holes;
  extract_contour_and_holes(alpha_shape, contour, holes);
  //! [Alpha Shapes]

  //! [Polyline Simplification]
  Cost cost;
  Stop stop(noise);
  CGAL::Polyline_simplification_2::simplify(contour, cost, stop);
  for (auto& hole : holes)
    CGAL::Polyline_simplification_2::simplify(hole, cost, stop);
  //! [Polyline Simplification]

  //! [Contour Regularization] 
  Contour_regularization_2 regularization(
    CGAL::parameters::
    min_angle(min_angle).
    min_length(min_length));

  regularization.estimate_principal_directions(contour);
  regularization.regularize(contour);
  for (auto& hole : holes)
    regularization.regularize(hole);
  //! [Contour Regularization]

  //! [LOD0 Export]
  LOD0 lod0;
  get_lod0(contour, holes, lod0);
  switch (export_type_0) {
    case Export_type::POLYGON_WITH_HOLES:
      export_lod0_as_polygon_with_holes(lod0, output);
      break;

    case Export_type::TRIANGLE_SOUP:
      export_lod0_as_triangle_soup(lod0, output);
      break;

    case Export_type::CITY_GML:
      export_lod0_as_city_gml(lod0, output);
      break;

    default:
      export_lod0_as_polygon_with_holes(lod0, output);
  }
  //! [LOD0 Export]

  //! [LOD1 Export]
  LOD1 lod1;
  if (use_average)
    get_lod1_avg(contour, holes, lod1);
  else 
    get_lod1_max(contour, holes, lod1);

  switch (export_type_1) {
    case Export_type::TRIANGLE_SOUP:
      export_lod1_as_triangle_soup(lod1, output);
      break;

    case Export_type::CITY_GML:
      export_lod1_as_city_gml(lod1, output);
      break;

    default:
      export_lod1_as_triangle_soup(lod1, output);
  }
  //! [LOD1 Export]

  return EXIT_SUCCESS;
}
