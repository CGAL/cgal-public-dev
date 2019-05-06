// CGAL includes.
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/IO/Color.h>
#include <CGAL/Point_set_3.h>
#include <CGAL/Point_set_3/IO.h>

// LOD includes.
#include <CGAL/Levels_of_detail.h>

// STL includes.
#include <iostream>

using Traits = CGAL::Exact_predicates_inexact_constructions_kernel;

using FT      = typename Traits::FT;
using Point_3 = typename Traits::Point_3;

using Points  = std::vector<Point_3>;
using Indices = std::vector< std::vector<std::size_t> >;
using Colors  = std::vector<CGAL::Color>;

using Point_set = Point_set_3<Point_3>;
using Point_map = typename Point_set::Point_map;
using Label_map = typename Point_set:: template Property_map<int>;

using Semantic_map   = Semantic_from_label_map<Label_map>;
using Visibility_map = Visibility_from_semantic_map<Semantic_map>;

using LOD = Levels_of_detail<
  Traits, 
  Point_set, 
  Point_map, 
  Semantic_map, 
  Visibility_map,
  CGAL::Tag_true>;

int main(int argc, char **argv) {

  // DS.
  Point_set point_set;
  Label_map label_map;

  // Load input data.
  std::cout << std::endl << "Input data: " << std::endl;
  std::ifstream file("data/lods.ply", std::ios_base::binary);
  file >> point_set; 
  file.close();
  std::cout << "File contains " << point_set.size() << " points" << std::endl;

  // Parameters.
  const FT ground_precision = FT(1);

  // Define a map from a user-defined label to the LOD semantic label.
  Semantic_map semantic_map(label_map, 0, 1, 2, 3);

  // Define a map for computing visibility.
  Visibility_map visibility_map(semantic_map);

  // Create LOD.
  LOD lod(
    point_set, 
    point_set.point_map(), 
    semantic_map,
    visibility_map);

  // Compute planar ground.
  Points pl_vertices; Indices pl_faces; Colors pl_fcolors;
  Polygon_inserter<Traits> pl_inserter(pl_faces, pl_fcolors);
  const auto success = lod.ground(
    std::back_inserter(pl_vertices),
    boost::make_function_output_iterator(pl_inserter),
    Reconstruction_type::PLANAR_GROUND,
    ground_precision);

  // Compute smooth ground.
  Points sm_vertices; Indices sm_faces; Colors sm_fcolors;
  Polygon_inserter<Traits> sm_inserter(sm_faces, sm_fcolors);
  const auto success = lod.ground(
    std::back_inserter(sm_vertices),
    boost::make_function_output_iterator(sm_inserter),
    Reconstruction_type::SMOOTH_GROUND,
    ground_precision);

  return EXIT_SUCCESS;
}
