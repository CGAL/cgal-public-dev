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
  const FT scale            = FT(4);
  const FT noise_level      = FT(2);
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

  // Compute all objects.
  lod.build(
    scale, noise_level);


  // Access LOD0.
  Points vertices_lod0; Indices_container faces_lod0; Colors fcolors_lod0;
  Polygon_inserter inserter_lod0(faces_lod0, fcolors_lod0);
  lod.lods(
    std::back_inserter(vertices_lod0),
    boost::make_function_output_iterator(inserter_lod0),
    CGAL::Levels_of_detail::Reconstruction_type::LOD0,
    ground_precision);

  // Access LOD1.
  Points vertices_lod1; Indices_container faces_lod1; Colors fcolors_lod1;
  Polygon_inserter inserter_lod1(faces_lod1, fcolors_lod1);
  lod.lods(
    std::back_inserter(vertices_lod1),
    boost::make_function_output_iterator(inserter_lod1),
    CGAL::Levels_of_detail::Reconstruction_type::LOD1,
    ground_precision);

  // Access LOD2.
  Points vertices_lod2; Indices_container faces_lod2; Colors fcolors_lod2;
  Polygon_inserter inserter_lod2(faces_lod2, fcolors_lod2);
  lod.lods(
    std::back_inserter(vertices_lod2),
    boost::make_function_output_iterator(inserter_lod2),
    CGAL::Levels_of_detail::Reconstruction_type::LOD2,
    ground_precision);


  std::cout << std::endl << std::endl;
  std::cout << "Number of lod0 faces: " << faces_lod0.size() << std::endl;
  std::cout << "Number of lod1 faces: " << faces_lod1.size() << std::endl;
  std::cout << "Number of lod2 faces: " << faces_lod2.size() << std::endl;

  return EXIT_SUCCESS;
}
