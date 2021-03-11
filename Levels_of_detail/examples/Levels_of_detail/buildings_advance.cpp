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
  const FT scale                                = FT(4);
  const FT noise_level                          = FT(2);
  const FT cluster_scale                        = FT(8);
  const std::size_t min_cluster_size            = 10;
  const FT alpha_shape_size_2                   = FT(2);
  const FT grid_cell_width_2                    = FT(1);
  const FT region_growing_scale_2               = FT(12);
  const FT region_growing_noise_level_2         = FT(2);
  const FT region_growing_angle_2               = FT(25);
  const FT region_growing_min_length_2          = FT(4);
  const FT kinetic_min_face_width_2             = FT(2);
  const std::size_t kinetic_max_intersections_2 = 2;
  const std::size_t min_faces_per_footprint     = 1;
  const FT visibility_scale_2                   = FT(4);
  const FT graphcut_beta_2                      = FT(1) / FT(10);
  const CGAL::Levels_of_detail::Extrusion_type extrusion_type =
  CGAL::Levels_of_detail::Extrusion_type::MAX;
  const FT region_growing_scale_3               = FT(2);
  const FT region_growing_noise_level_3         = FT(2);
  const FT region_growing_angle_3               = FT(25);
  const FT region_growing_min_area_3            = FT(8);
  const FT region_growing_distance_to_line_3    = FT(1);
  const std::size_t kinetic_max_intersections_3 = 2;
  const FT visibility_scale_3                   = FT(4);
  const FT graphcut_beta_3                      = FT(1) / FT(10);

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
  Point_inserter inserter_ibb(points);
  lod.points(
    boost::make_function_output_iterator(inserter_ibb),
    CGAL::Levels_of_detail::Intermediate_step::INPUT_BUILDING_BOUNDARY_POINTS);
  std::cout << "Number of building boundary points: " << points.size()
  << std::endl;

  points.clear();
  Point_inserter inserter_ibi(points);
  lod.points(
    boost::make_function_output_iterator(inserter_ibi),
    CGAL::Levels_of_detail::Intermediate_step::INPUT_BUILDING_INTERIOR_POINTS);
  std::cout << "Number of building interior points: " << points.size()
  << std::endl << std::endl;


  // Initialize buildings - STEP 1.
  lod.initialize_buildings(
    scale, noise_level, cluster_scale, min_cluster_size);

  // Access intermediate steps.
  points.clear();
  Point_inserter inserter_bcl(points);
  lod.points(
  boost::make_function_output_iterator(inserter_bcl),
  CGAL::Levels_of_detail::Intermediate_step::BUILDING_CLUSTERS);
  std::cout << "Number of cluster points: " << points.size() << std::endl;


  // Detect building boundaries - STEP 2.
  lod.detect_building_boundaries(
    alpha_shape_size_2, grid_cell_width_2,
    region_growing_scale_2, region_growing_noise_level_2,
    region_growing_angle_2, region_growing_min_length_2);

  // Access intermediate steps.
  points.clear();
  Point_inserter inserter_bbp(points);
  lod.points(
  boost::make_function_output_iterator(inserter_bbp),
  CGAL::Levels_of_detail::Intermediate_step::BUILDING_BOUNDARY_POINTS);
  std::cout << "Number of boundary points: " << points.size() << std::endl;

  points.clear();
  Point_inserter inserter_bwp(points);
  lod.points(
  boost::make_function_output_iterator(inserter_bwp),
  CGAL::Levels_of_detail::Intermediate_step::BUILDING_WALL_POINTS);
  std::cout << "Number of wall points: " << points.size() << std::endl;

  Points_container segments;
  Polyline_inserter inserter_bab(segments);
  lod.polylines(
    boost::make_function_output_iterator(inserter_bab),
    CGAL::Levels_of_detail::Intermediate_step::BUILDING_APPROXIMATE_BOUNDARIES);
  std::cout << "Number of approximate boundary edges: " << segments.size() << std::endl;


  // Compute building footprints - STEP 3.
  lod.compute_building_footprints(
    kinetic_min_face_width_2, kinetic_max_intersections_2,
    min_faces_per_footprint, visibility_scale_2, graphcut_beta_2);

  // Access intermediate steps.
  Points vertices; Indices_container faces; Colors fcolors;
  Polygon_inserter inserter_bp2(faces, fcolors);
  lod.mesh(
    std::back_inserter(vertices),
    boost::make_function_output_iterator(inserter_bp2),
    CGAL::Levels_of_detail::Intermediate_step::BUILDING_PARTITIONING_2);
  std::cout << "Number of partitioning2 faces: " << faces.size() << std::endl;

  points.clear();
  Point_inserter inserter_bpo(points);
  lod.points(
  boost::make_function_output_iterator(inserter_bpo),
  CGAL::Levels_of_detail::Intermediate_step::BUILDING_POINTS);
  std::cout << "Number of building points: " << points.size() << std::endl;

  segments.clear();
  Polyline_inserter inserter_bbo(segments);
  lod.polylines(
    boost::make_function_output_iterator(inserter_bbo),
    CGAL::Levels_of_detail::Intermediate_step::BUILDING_BOUNDARIES);
  std::cout << "Number of boundary edges: " << segments.size() << std::endl;

  vertices.clear(); faces.clear(); fcolors.clear();
  Polygon_inserter inserter_bfo(faces, fcolors);
  lod.mesh(
    std::back_inserter(vertices),
    boost::make_function_output_iterator(inserter_bfo),
    CGAL::Levels_of_detail::Intermediate_step::BUILDING_FOOTPRINTS);
  std::cout << "Number of footprint faces: " << faces.size() << std::endl;


  // Extrude building footprints - STEP 4.
  lod.extrude_building_footprints(
    extrusion_type);

  // Access intermediate steps.
  vertices.clear(); faces.clear(); fcolors.clear();
  Polygon_inserter inserter_ebb(faces, fcolors);
  lod.mesh(
    std::back_inserter(vertices),
    boost::make_function_output_iterator(inserter_ebb),
    CGAL::Levels_of_detail::Intermediate_step::EXTRUDED_BUILDING_BOUNDARIES);
  std::cout << "Number of extruded boundary faces: "
  << faces.size() << std::endl;

  vertices.clear(); faces.clear(); fcolors.clear();
  Polygon_inserter inserter_ebf(faces, fcolors);
  lod.mesh(
    std::back_inserter(vertices),
    boost::make_function_output_iterator(inserter_ebf),
    CGAL::Levels_of_detail::Intermediate_step::EXTRUDED_BUILDING_FOOTPRINTS);
  std::cout << "Number of extruded footprint faces: "
  << faces.size() << std::endl;


  // Detect building roofs - STEP 5.
  lod.detect_building_roofs(
    region_growing_scale_3, region_growing_noise_level_3,
    region_growing_angle_3, region_growing_min_area_3,
    region_growing_distance_to_line_3);

  // Access intermediate steps.
  points.clear();
  Point_inserter inserter_brp(points);
  lod.points(
  boost::make_function_output_iterator(inserter_brp),
  CGAL::Levels_of_detail::Intermediate_step::BUILDING_ROOF_POINTS);
  std::cout << "Number of roof points: " << points.size() << std::endl;

  vertices.clear(); faces.clear(); fcolors.clear();
  Polygon_inserter inserter_abb1(faces, fcolors);
  lod.mesh(
    std::back_inserter(vertices),
    boost::make_function_output_iterator(inserter_abb1),
    CGAL::Levels_of_detail::Intermediate_step::APPROXIMATE_BUILDING_BOUNDS);
  std::cout << "Number of approximate bound faces 1: "
  << faces.size() << std::endl;


  // Compute building roofs - STEP 6.
  lod.compute_building_roofs(
    kinetic_max_intersections_3,
    visibility_scale_3, graphcut_beta_3);

  // Access intermediate steps.
  vertices.clear(); faces.clear(); fcolors.clear();
  Polygon_inserter inserter_abb2(faces, fcolors);
  lod.mesh(
    std::back_inserter(vertices),
    boost::make_function_output_iterator(inserter_abb2),
    CGAL::Levels_of_detail::Intermediate_step::APPROXIMATE_BUILDING_BOUNDS);
  std::cout << "Number of approximate bound faces 2: "
  << faces.size() << std::endl;

  vertices.clear(); faces.clear(); fcolors.clear();
  Polygon_inserter inserter_bp3(faces, fcolors);
  lod.mesh(
    std::back_inserter(vertices),
    boost::make_function_output_iterator(inserter_bp3),
    CGAL::Levels_of_detail::Intermediate_step::BUILDING_PARTITIONING_3);
  std::cout << "Number of partitioning3 faces: " << faces.size() << std::endl;

  vertices.clear(); faces.clear(); fcolors.clear();
  Polygon_inserter inserter_bwa(faces, fcolors);
  lod.mesh(
    std::back_inserter(vertices),
    boost::make_function_output_iterator(inserter_bwa),
    CGAL::Levels_of_detail::Intermediate_step::BUILDING_WALLS);
  std::cout << "Number of wall faces: "
  << faces.size() << std::endl;

  vertices.clear(); faces.clear(); fcolors.clear();
  Polygon_inserter inserter_bro(faces, fcolors);
  lod.mesh(
    std::back_inserter(vertices),
    boost::make_function_output_iterator(inserter_bro),
    CGAL::Levels_of_detail::Intermediate_step::BUILDING_ROOFS);
  std::cout << "Number of roof faces: "
  << faces.size() << std::endl;


  // Access BUILDINGS0.
  std::cout << std::endl << std::endl;
  vertices.clear(); faces.clear(); fcolors.clear();
  Polygon_inserter inserter_bu0(faces, fcolors);
  lod.buildings(
    std::back_inserter(vertices),
    boost::make_function_output_iterator(inserter_bu0),
    CGAL::Levels_of_detail::Reconstruction_type::BUILDINGS0);
  std::cout << "Number of buildings0 faces: " << faces.size() << std::endl;

  // Access BUILDINGS1.
  vertices.clear(); faces.clear(); fcolors.clear();
  Polygon_inserter inserter_bu1(faces, fcolors);
  lod.buildings(
    std::back_inserter(vertices),
    boost::make_function_output_iterator(inserter_bu1),
    CGAL::Levels_of_detail::Reconstruction_type::BUILDINGS1);
  std::cout << "Number of buildings1 faces: " << faces.size() << std::endl;

  // Access BUILDINGS2.
  vertices.clear(); faces.clear(); fcolors.clear();
  Polygon_inserter inserter_bu2(faces, fcolors);
  lod.buildings(
    std::back_inserter(vertices),
    boost::make_function_output_iterator(inserter_bu2),
    CGAL::Levels_of_detail::Reconstruction_type::BUILDINGS2);
  std::cout << "Number of buildings2 faces: " << faces.size() << std::endl;

  return EXIT_SUCCESS;
}
