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
  const FT scale       = FT(4);
  const FT noise_level = FT(2);

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

  // Compute trees.
  lod.build_trees(
    scale, noise_level);


  // Access TREES0.
  std::cout << std::endl << std::endl;
  Points vertices; Indices_container faces; Colors fcolors;
  Polygon_inserter inserter_tr0(faces, fcolors);
  lod.trees(
    std::back_inserter(vertices),
    boost::make_function_output_iterator(inserter_tr0),
    CGAL::Levels_of_detail::Reconstruction_type::TREES0);
  std::cout << "Number of trees0 faces: " << faces.size() << std::endl;

  // Access TREES1.
  vertices.clear(); faces.clear(); fcolors.clear();
  Polygon_inserter inserter_tr1(faces, fcolors);
  lod.trees(
    std::back_inserter(vertices),
    boost::make_function_output_iterator(inserter_tr1),
    CGAL::Levels_of_detail::Reconstruction_type::TREES1);
  std::cout << "Number of trees1 faces: " << faces.size() << std::endl;

  // Access TREES2.
  vertices.clear(); faces.clear(); fcolors.clear();
  Polygon_inserter inserter_tr2(faces, fcolors);
  lod.trees(
    std::back_inserter(vertices),
    boost::make_function_output_iterator(inserter_tr2),
    CGAL::Levels_of_detail::Reconstruction_type::TREES2);
  std::cout << "Number of trees2 faces: " << faces.size() << std::endl;

  return EXIT_SUCCESS;
}
