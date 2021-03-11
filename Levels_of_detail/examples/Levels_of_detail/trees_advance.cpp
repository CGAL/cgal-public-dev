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
using Points_container  = std::vector<Points>;
using Indices           = std::vector<std::size_t>;
using Indices_container = std::vector<Indices>;
using Colors            = std::vector<CGAL::Color>;

using Point_set = CGAL::Point_set_3<Point_3>;
using Point_map = typename Point_set::Point_map;
using Label_map = typename Point_set:: template Property_map<int>;

using Semantic_map   = CGAL::Levels_of_detail::Semantic_from_label_map<Label_map>;
using Visibility_map = CGAL::Levels_of_detail::Visibility_from_semantic_map<Semantic_map>;

using Point_inserter    = CGAL::Levels_of_detail::Point_inserter<Traits>;
using Polygon_inserter  = CGAL::Levels_of_detail::Polygon_inserter<Traits>;
using Polyline_inserter = CGAL::Levels_of_detail::Polyline_inserter<Traits>;

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
  const FT scale                            = FT(4);
  const FT noise_level                      = FT(2);
  const FT cluster_scale                    = FT(8);
  const std::size_t min_cluster_size        = 10;
  const FT grid_cell_width_2                = FT(1) / FT(2);
  const FT min_height                       = FT(3);
  const FT min_radius_2                     = FT(2);
  const std::size_t min_faces_per_footprint = 12;
  const CGAL::Levels_of_detail::Extrusion_type extrusion_type =
  CGAL::Levels_of_detail::Extrusion_type::MAX;

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
  Point_inserter inserter_ivp(points);
  lod.points(
    boost::make_function_output_iterator(inserter_ivp),
    CGAL::Levels_of_detail::Intermediate_step::INPUT_VEGETATION_POINTS);
  std::cout << "Number of vegetation points: " << points.size()
  << std::endl << std::endl;


  // Initialize trees - STEP 1.
  lod.initialize_trees(
    scale, noise_level, cluster_scale, min_cluster_size);

  // Access intermediate steps.
  points.clear();
  Point_inserter inserter_tcl(points);
  lod.points(
  boost::make_function_output_iterator(inserter_tcl),
  CGAL::Levels_of_detail::Intermediate_step::TREE_CLUSTERS);
  std::cout << "Number of cluster points: " << points.size() << std::endl;


  // Compute tree footprints - STEP 2.
  lod.compute_tree_footprints(
    grid_cell_width_2, min_height, min_radius_2, min_faces_per_footprint);

  // Access intermediate steps.
  points.clear();
  Point_inserter inserter_tpo(points);
  lod.points(
  boost::make_function_output_iterator(inserter_tpo),
  CGAL::Levels_of_detail::Intermediate_step::TREE_POINTS);
  std::cout << "Number of tree points: " << points.size() << std::endl;

  Points_container segments;
  Polyline_inserter inserter_tbo(segments);
  lod.polylines(
    boost::make_function_output_iterator(inserter_tbo),
    CGAL::Levels_of_detail::Intermediate_step::TREE_BOUNDARIES);
  std::cout << "Number of boundary edges: " << segments.size() << std::endl;

  Points vertices; Indices_container faces; Colors fcolors;
  Polygon_inserter inserter_tfo(faces, fcolors);
  lod.mesh(
    std::back_inserter(vertices),
    boost::make_function_output_iterator(inserter_tfo),
    CGAL::Levels_of_detail::Intermediate_step::TREE_FOOTPRINTS);
  std::cout << "Number of footprint faces: " << faces.size() << std::endl;


  // Extrude tree footprints - STEP 3.
  lod.extrude_tree_footprints(
    extrusion_type);

  // Access intermediate steps.
  vertices.clear(); faces.clear(); fcolors.clear();
  Polygon_inserter inserter_etb(faces, fcolors);
  lod.mesh(
    std::back_inserter(vertices),
    boost::make_function_output_iterator(inserter_etb),
    CGAL::Levels_of_detail::Intermediate_step::EXTRUDED_TREE_BOUNDARIES);
  std::cout << "Number of extruded boundary faces: "
  << faces.size() << std::endl;

  vertices.clear(); faces.clear(); fcolors.clear();
  Polygon_inserter inserter_etf(faces, fcolors);
  lod.mesh(
    std::back_inserter(vertices),
    boost::make_function_output_iterator(inserter_etf),
    CGAL::Levels_of_detail::Intermediate_step::EXTRUDED_TREE_FOOTPRINTS);
  std::cout << "Number of extruded footprint faces: "
  << faces.size() << std::endl;


  // Compute tree crowns - STEP 4.
  lod.compute_tree_crowns();

  // Access intermediate steps.
  vertices.clear(); faces.clear(); fcolors.clear();
  Polygon_inserter inserter_ttr(faces, fcolors);
  lod.mesh(
    std::back_inserter(vertices),
    boost::make_function_output_iterator(inserter_ttr),
    CGAL::Levels_of_detail::Intermediate_step::TREE_TRUNKS);
  std::cout << "Number of trunk faces: "
  << faces.size() << std::endl;

  vertices.clear(); faces.clear(); fcolors.clear();
  Polygon_inserter inserter_tcr(faces, fcolors);
  lod.mesh(
    std::back_inserter(vertices),
    boost::make_function_output_iterator(inserter_tcr),
    CGAL::Levels_of_detail::Intermediate_step::TREE_CROWNS);
  std::cout << "Number of crown faces: "
  << faces.size() << std::endl;


  // Access TREES0.
  std::cout << std::endl << std::endl;
  vertices.clear(); faces.clear(); fcolors.clear();
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
