// STL includes.
#include <vector>
#include <fstream>
#include <iostream>

// CGAL includes.
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/IO/Color.h>
#include <CGAL/Point_set_3.h>
#include <CGAL/Point_set_3/IO.h>

// LOD includes.
#include <CGAL/Levels_of_detail.h>

// Internal includes.
#include "../../test/Levels_of_detail/include/Utilities.h"

using Traits = CGAL::Exact_predicates_inexact_constructions_kernel;

using FT      = typename Traits::FT;
using Point_3 = typename Traits::Point_3;

using Points            = std::vector<Point_3>;
using Indices           = std::vector<std::size_t>;
using Indices_container = std::vector<Indices>;
using Colors            = std::vector<CGAL::Color>;

using Point_set = CGAL::Point_set_3<Point_3>;
using Point_map = typename Point_set::Point_map;
using Label_map = typename Point_set:: template Property_map<int>;

using Semantic_map   = CGAL::Levels_of_detail::Semantic_from_label_map<Label_map>;
using Visibility_map = CGAL::Levels_of_detail::Visibility_from_semantic_map<Semantic_map>;

using Point_inserter   = CGAL::Levels_of_detail::Point_inserter<Traits>;
using Polygon_inserter = CGAL::Levels_of_detail::Polygon_inserter<Traits>;

using LOD = CGAL::Levels_of_detail::Levels_of_detail<
  Traits,
  Point_set,
  Point_map,
  Semantic_map,
  Visibility_map,
  CGAL::Tag_true>;

int main(int argc, char **argv) {

  // Load input data.
  Point_set point_set;
  std::cout << std::endl << "Input data: " << std::endl;
  std::ifstream file(argc > 1 ? argv[1] : "data/lods.ply",
  std::ios_base::binary);
  file >> point_set;
  file.close();
  std::cout << "File contains " << point_set.size() << " points" << std::endl;

  // Parameters.
  const FT ground_precision = FT(4);

  // Define a map from a user-defined label to the LOD semantic label.
  Label_map label_map = point_set. template property_map<int>("label").first;
  Semantic_map semantic_map(label_map, "0", "1", "2", "3");

  // Define a map for computing visibility.
  Visibility_map visibility_map(semantic_map);


  // Initialize LOD.
  LOD lod(
    point_set,
    point_set.point_map(),
    semantic_map,
    visibility_map);

  // Access intermediate steps.
  Point_set points;
  Point_inserter inserter_igp(points);
  lod.points(
    boost::make_function_output_iterator(inserter_igp),
    CGAL::Levels_of_detail::Intermediate_step::INPUT_GROUND_POINTS);
  std::cout << "Number of ground points: " << points.size()
  << std::endl << std::endl;


  // Compute planar ground.
  Points vertices_pl; Indices_container faces_pl; Colors fcolors_pl;
  Polygon_inserter inserter_pl(faces_pl, fcolors_pl);
  lod.ground(
    std::back_inserter(vertices_pl),
    boost::make_function_output_iterator(inserter_pl),
    CGAL::Levels_of_detail::Reconstruction_type::PLANAR_GROUND,
    ground_precision);

  // Compute smooth ground.
  Points vertices_sm; Indices_container faces_sm; Colors fcolors_sm;
  Polygon_inserter inserter_sm(faces_sm, fcolors_sm);
  lod.ground(
    std::back_inserter(vertices_sm),
    boost::make_function_output_iterator(inserter_sm),
    CGAL::Levels_of_detail::Reconstruction_type::SMOOTH_GROUND,
    ground_precision);


  std::cout << std::endl << std::endl;
  std::cout << "Number of planar ground faces: " << faces_pl.size() << std::endl;
  std::cout << "Number of smooth ground faces: " << faces_sm.size() << std::endl;

  return EXIT_SUCCESS;
}
