#ifndef CGAL_LOD_WRAPPER_H
#define CGAL_LOD_WRAPPER_H

#if defined(WIN32) || defined(_WIN32)
#define _SR_ "\\"
#else
#define _SR_ "/"
#endif

// STL includes.
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>

// Boost includes.
#include <boost/function_output_iterator.hpp>

// CGAL includes.
#include <CGAL/Timer.h>
#include <CGAL/IO/Color.h>
#include <CGAL/Point_set_3.h>
#include <CGAL/Point_set_3/IO.h>

// LOD includes.
#include <CGAL/Levels_of_detail.h>

// Local includes.
#include "Saver.h"
#include "Utilities.h"
#include "Terminal_parser.h"

namespace CGAL {
namespace Levels_of_detail {

  template<typename GeomTraits>
  class Wrapper {

  public:
    using Traits = GeomTraits;

    using FT = typename Traits::FT;
    using Point_3 = typename Traits::Point_3;
    using Vector_3 = typename Traits::Vector_3;

    using Saver = Saver<Traits>;
    using Parameters = internal::Parameters<FT>;
    using Terminal_parser = Terminal_parser<FT>;
    using Point_set = Point_set_3<Point_3>;

    using Points = std::vector<Point_3>;
    using Points_container = std::vector<Points>;
    using Indices = std::vector<std::size_t>;
    using Indices_container = std::vector<Indices>;
    using Colors = std::vector<CGAL::Color>;

    using Point_map = typename Point_set::Point_map;
    using Label_map = typename Point_set:: template Property_map<int>;

    using Semantic_map = Semantic_from_label_map<Label_map>;
    using Visibility_map = Visibility_from_semantic_map<Semantic_map>;

    using LOD = Levels_of_detail<
    Traits,
    Point_set,
    Point_map,
    Semantic_map,
    Visibility_map,
    CGAL::Tag_true>;

    Wrapper(
      int argc,
      char **argv,
      const std::string path_to_save) :
    m_terminal_parser(argc, argv, path_to_save),
    m_path(path_to_save),
    m_path_gr(m_path + "ground" + std::string(_SR_)),
    m_path_tr(m_path + "trees" + std::string(_SR_)),
    m_path_bu(m_path + "buildings" + std::string(_SR_)),
    m_make_ground(false),
    m_make_trees(false),
    m_make_buildings(false),
    m_make_lod01(false),
    m_make_lod2(false),
    m_make_all(false)
    { }

    void execute() {
      parse_terminal();
      load_input_data();
      execute_pipeline();
    }

  private:
    Saver m_saver;
    Parameters m_parameters;
    Terminal_parser m_terminal_parser;
    std::string m_path, m_path_gr, m_path_tr, m_path_bu;
    Point_set m_point_set;
    Label_map m_label_map;

    bool m_make_ground;
    bool m_make_trees;
    bool m_make_buildings;
    bool m_make_lod01;
    bool m_make_lod2;
    bool m_make_all;

    void parse_terminal() {
      // Set all parameters that can be loaded from the terminal.
      // add_str_parameter  - adds a string-type parameter
      // add_val_parameter  - adds a scalar-type parameter
      // add_bool_parameter - adds a boolean parameter

      std::cout << "Input parameters: " << std::endl;

      // Required parameters.
      m_terminal_parser.add_str_parameter("-data", m_parameters.data);

      // Label indices.
      m_terminal_parser.add_str_parameter("-gi", m_parameters.gi);
      m_terminal_parser.add_str_parameter("-bi", m_parameters.bi);
      m_terminal_parser.add_str_parameter("-ii", m_parameters.ii);
      m_terminal_parser.add_str_parameter("-vi", m_parameters.vi);

      // Main parameters.
      m_terminal_parser.add_val_parameter("-scale", m_parameters.scale);
      m_terminal_parser.add_val_parameter("-noise", m_parameters.noise_level);

      m_terminal_parser.add_bool_parameter("-lidar", m_parameters.lidar);

      m_terminal_parser.add_bool_parameter("-make_ground", m_make_ground);
      m_terminal_parser.add_bool_parameter("-make_trees", m_make_trees);
      m_terminal_parser.add_bool_parameter("-make_buildings", m_make_buildings);
      m_terminal_parser.add_bool_parameter("-make_lod01", m_make_lod01);
      m_terminal_parser.add_bool_parameter("-make_lod2", m_make_lod2);
      m_terminal_parser.add_bool_parameter("-make_all", m_make_all);


      // Update.
      m_parameters.update_dependent();


      // Clustering buildings.
      m_terminal_parser.add_val_parameter("-bu_clust", m_parameters.buildings.cluster_scale);
      m_terminal_parser.add_val_parameter("-bu_min", m_parameters.buildings.min_cluster_size);

      // Detecting building boundaries.
      m_terminal_parser.add_val_parameter("-alpha_2", m_parameters.buildings.alpha_shape_size_2);
      m_terminal_parser.add_val_parameter("-bu_cell_2", m_parameters.buildings.grid_cell_width_2);
      m_terminal_parser.add_val_parameter("-im_beta_2", m_parameters.buildings.imagecut_beta_2);
      m_terminal_parser.add_val_parameter("-im_noise_2", m_parameters.buildings.image_noise_2);

      m_terminal_parser.add_val_parameter("-rg_scale_2", m_parameters.buildings.region_growing_scale_2);
      m_terminal_parser.add_val_parameter("-rg_noise_2", m_parameters.buildings.region_growing_noise_level_2);
      m_terminal_parser.add_val_parameter("-rg_angle_2", m_parameters.buildings.region_growing_angle_2);
      m_terminal_parser.add_val_parameter("-rg_length_2", m_parameters.buildings.region_growing_min_length_2);

      m_terminal_parser.add_val_parameter("-reg_length_2", m_parameters.buildings.regularization_min_length_2);
      m_terminal_parser.add_val_parameter("-reg_angle_2", m_parameters.buildings.regularization_angle_bound_2);
      m_terminal_parser.add_val_parameter("-reg_ordinate_2", m_parameters.buildings.regularization_ordinate_bound_2);

      // Computing building footprints.
      m_terminal_parser.add_val_parameter("-bu_faces", m_parameters.buildings.min_faces_per_footprint);
      m_terminal_parser.add_val_parameter("-kn_width_2", m_parameters.buildings.kinetic_min_face_width_2);
      m_terminal_parser.add_val_parameter("-kn_inter_2", m_parameters.buildings.kinetic_max_intersections_2);
      m_terminal_parser.add_val_parameter("-vis_scale_2", m_parameters.buildings.visibility_scale_2);
      m_terminal_parser.add_val_parameter("-gc_beta_2", m_parameters.buildings.graphcut_beta_2);

      // Detecting building roofs.
      m_terminal_parser.add_val_parameter("-rg_scale_3", m_parameters.buildings.region_growing_scale_3);
      m_terminal_parser.add_val_parameter("-rg_noise_3", m_parameters.buildings.region_growing_noise_level_3);
      m_terminal_parser.add_val_parameter("-rg_angle_3", m_parameters.buildings.region_growing_angle_3);
      m_terminal_parser.add_val_parameter("-rg_area_3", m_parameters.buildings.region_growing_min_area_3);
      m_terminal_parser.add_val_parameter("-rg_dist_3", m_parameters.buildings.region_growing_distance_to_line_3);
      m_terminal_parser.add_val_parameter("-max_height", m_parameters.buildings.max_height_difference);

      // Computing building roofs.
      m_terminal_parser.add_val_parameter("-kn_inter_3", m_parameters.buildings.kinetic_max_intersections_3);
      m_terminal_parser.add_val_parameter("-vis_scale_3", m_parameters.buildings.visibility_scale_3);
      m_terminal_parser.add_val_parameter("-gc_beta_3", m_parameters.buildings.graphcut_beta_3);


      // Clustering trees.
      m_terminal_parser.add_val_parameter("-tr_clust", m_parameters.trees.cluster_scale);
      m_terminal_parser.add_val_parameter("-tr_min", m_parameters.trees.min_cluster_size);

      // Computing tree footprints.
      m_terminal_parser.add_val_parameter("-tr_cell_2", m_parameters.trees.grid_cell_width_2);
      m_terminal_parser.add_val_parameter("-tr_height", m_parameters.trees.min_height);
      m_terminal_parser.add_val_parameter("-tr_radius_2", m_parameters.trees.min_radius_2);
      m_terminal_parser.add_val_parameter("-tr_faces", m_parameters.trees.min_faces_per_footprint);


      // Smooth ground.
      m_terminal_parser.add_val_parameter("-gr_prec", m_parameters.ground.precision);
    }

    void load_input_data() {
      std::cout << std::endl << "Input data: " << std::endl;
      std::ifstream file(m_parameters.data.c_str(), std::ios_base::binary);
      file >> m_point_set;
      file.close();
      std::cout << "File contains " << m_point_set.size()
      << " points" << std::endl;

      if (are_label_data_defined()) {
        std::cout <<
          "Label data are defined!"
        << std::endl << std::endl;
        m_label_map = m_point_set. template property_map<int>("label").first;
      } else {
        std::cerr <<
          "Label data are not defined!"
        << std::endl << std::endl;
        exit(EXIT_FAILURE);
      }
    }

    bool are_label_data_defined() const {
      return m_point_set. template property_map<int>("label").second;
    }

    void execute_pipeline() {

      if (m_make_all) {
        m_make_ground = true;
        m_make_trees = true;
        m_make_buildings = true;
        m_make_lod01 = true;
        m_make_lod2 = true;
      }

      // Define a map from a user-defined label to the LOD semantic label.
      Semantic_map semantic_map(m_label_map,
      m_parameters.gi,
      m_parameters.bi,
      m_parameters.ii,
      m_parameters.vi);

      // Define a map for computing visibility.
      Visibility_map visibility_map(semantic_map);

      /*
      recenter_point_cloud();
      exit(EXIT_SUCCESS); */

      // Create LOD.
      LOD lod(
        m_point_set,
        m_point_set.point_map(),
        semantic_map,
        m_parameters,
        visibility_map,
        Input_type::POINT_SET);

      /*
      std::cout << std::endl << "COMPLETE:" << std::endl;
      lod.build();
      lod.build_trees();
      lod.build_buildings(); */

      std::cout << std::endl << "STEPS:" << std::endl;


      // Ground.
      if (m_make_ground) {
        save_ground_input(lod);

        if (m_make_lod01 || m_make_lod2) {
          save_ground(lod,
          Reconstruction_type::PLANAR_GROUND, m_parameters.ground.precision,
          m_path + "planar_ground");
          save_ground(lod,
          Reconstruction_type::SMOOTH_GROUND, m_parameters.ground.precision,
          m_path + "smooth_ground");
        }
      }


      // Trees.
      if (m_make_trees) {
        save_trees_input(lod);

        if (m_make_lod01) {
          lod.initialize_trees();
          save_tree_clusters(lod);

          lod.compute_tree_footprints();
          save_trees_before_extrusion(lod);

          lod.extrude_tree_footprints();
          save_trees_after_extrusion(lod);

          save_trees(lod, Reconstruction_type::TREES0, m_path + "trees0");
          save_trees(lod, Reconstruction_type::TREES1, m_path + "trees1");
        }

        if (m_make_lod2) {
          lod.compute_tree_crowns();
          save_trees_with_crowns(lod);

          save_trees(lod, Reconstruction_type::TREES2, m_path + "trees2");
        }
      }


      // Buildings.
      if (m_make_buildings) {
        save_buildings_input(lod);

        if (m_make_lod01) {
          lod.initialize_buildings();
          save_building_clusters(lod);

          lod.detect_building_boundaries();
          save_buildings_before_extrusion1(lod);

          lod.compute_building_footprints();
          save_buildings_before_extrusion2(lod);

          lod.extrude_building_footprints();
          save_buildings_after_extrusion(lod);

          save_buildings(lod, Reconstruction_type::BUILDINGS0, m_path + "buildings0");
          save_buildings(lod, Reconstruction_type::BUILDINGS1, m_path + "buildings1");
        }

        if (m_make_lod2) {
          lod.detect_building_roofs();
          save_roofs_before_extraction(lod);

          lod.compute_building_roofs();
          save_roofs_after_extraction(lod);

          save_buildings(lod, Reconstruction_type::BUILDINGS2, m_path + "buildings2");
        }
      }


      // LODs.
      if (m_make_lod01) {
        save_lod(lod,
        Reconstruction_type::LOD0, m_parameters.ground.precision,
        m_path + "LOD0");
        save_lod(lod,
        Reconstruction_type::LOD1, m_parameters.ground.precision,
        m_path + "LOD1");
      }

      if (m_make_lod2) {
        save_lod(lod,
        Reconstruction_type::LOD2, m_parameters.ground.precision,
        m_path + "LOD2");
      }
    }

    // Results.
    void save_ground(
      const LOD& lod,
      const Reconstruction_type ground_type,
      const FT ground_precision,
      const std::string path) {

      Points vertices; Indices_container faces; Colors fcolors;
      Polygon_inserter<Traits> inserter(faces, fcolors);

      const auto success = lod.ground(
        std::back_inserter(vertices),
        boost::make_function_output_iterator(inserter),
        ground_type,
        ground_precision);

      if (success)
        m_saver.export_polygon_soup(vertices, faces, fcolors, path);
    }

    void save_trees(
      const LOD& lod,
      const Reconstruction_type lod_type,
      const std::string path) {

      Points vertices; Indices_container faces; Colors fcolors;
      Polygon_inserter<Traits> inserter(faces, fcolors);

      const auto success = lod.trees(
        std::back_inserter(vertices),
        boost::make_function_output_iterator(inserter),
        lod_type);

      if (success)
        m_saver.export_polygon_soup(vertices, faces, fcolors, path);
    }

    void save_buildings(
      const LOD& lod,
      const Reconstruction_type lod_type,
      const std::string path) {

      Points vertices; Indices_container faces; Colors fcolors;
      Polygon_inserter<Traits> inserter(faces, fcolors);

      const auto success = lod.buildings(
        std::back_inserter(vertices),
        boost::make_function_output_iterator(inserter),
        lod_type);

      if (success)
        m_saver.export_polygon_soup(vertices, faces, fcolors, path);
    }

    void save_lod(
      const LOD& lod,
      Reconstruction_type lod_type,
      const FT ground_precision,
      const std::string path) {

      Points vertices; Indices_container faces; Colors fcolors;
      Polygon_inserter<Traits> inserter(faces, fcolors);

      const auto success = lod.lods(
        std::back_inserter(vertices),
        boost::make_function_output_iterator(inserter),
        lod_type,
        ground_precision);

      if (success)
        m_saver.export_polygon_soup(vertices, faces, fcolors, path);
    }

    // Inermediate ground.
    void save_ground_input(const LOD& lod) {
      save_points(lod,
      Intermediate_step::INPUT_GROUND_POINTS,
      m_path_gr + "0_ground_points");
    }

    // Inermediate trees.
    void save_trees_input(const LOD& lod) {
      save_points(lod,
      Intermediate_step::INPUT_VEGETATION_POINTS,
      m_path_tr + "0_vegetation_points");
    }

    void save_tree_clusters(const LOD& lod) {
      save_points(lod, Intermediate_step::TREE_CLUSTERS,
      m_path_tr + "trees_1_clusters");
    }

    void save_trees_before_extrusion(const LOD& lod) {
      save_points(lod, Intermediate_step::TREE_POINTS,
      m_path_tr + "trees_2_points");
      save_polylines(lod, Intermediate_step::TREE_BOUNDARIES,
      m_path_tr + "trees_3_boundaries");
      save_mesh(lod, Intermediate_step::TREE_FOOTPRINTS,
      m_path_tr + "trees_4_footprints");
    }

    void save_trees_after_extrusion(const LOD& lod) {
      save_mesh(lod, Intermediate_step::EXTRUDED_TREE_BOUNDARIES,
      m_path_tr + "trees_5_extruded_boundaries");
      save_mesh(lod, Intermediate_step::EXTRUDED_TREE_FOOTPRINTS,
      m_path_tr + "trees_6_extruded_footprints");
    }

    void save_trees_with_crowns(const LOD& lod) {
      save_mesh(lod, Intermediate_step::TREE_TRUNKS,
      m_path_tr + "trees_7_trunks");
      save_mesh(lod, Intermediate_step::TREE_CROWNS,
      m_path_tr + "trees_8_crowns");
    }

    // Inermediate buildings.
    void save_buildings_input(const LOD& lod) {
      save_points(lod,
      Intermediate_step::INPUT_BUILDING_BOUNDARY_POINTS,
      m_path_bu + "01_building_boundary_points");
      save_points(lod,
      Intermediate_step::INPUT_BUILDING_INTERIOR_POINTS,
      m_path_bu + "02_building_interior_points");
    }

    void save_building_clusters(const LOD& lod) {
      save_points(lod, Intermediate_step::BUILDING_CLUSTERS,
      m_path_bu + "buildings_1_clusters");
    }

    void save_buildings_before_extrusion1(const LOD& lod) {
      save_points(lod, Intermediate_step::BUILDING_BOUNDARY_POINTS,
      m_path_bu + "buildings_2_boundary_points");
      /*
      save_points(lod, Intermediate_step::BUILDING_WALL_POINTS,
      m_path_bu + "buildings_2_wall_points"); */
      save_polylines(lod, Intermediate_step::BUILDING_APPROXIMATE_BOUNDARIES,
      m_path_bu + "buildings_3_approximate_boundaries");
    }

    void save_buildings_before_extrusion2(const LOD& lod) {
      save_mesh(lod, Intermediate_step::BUILDING_PARTITIONING_2,
      m_path_bu + "buildings_4_partitioning_2");
      save_polylines(lod, Intermediate_step::BUILDING_BOUNDARIES,
      m_path_bu + "buildings_5_boundaries");
      save_mesh(lod, Intermediate_step::BUILDING_FOOTPRINTS,
      m_path_bu + "buildings_6_footprints");
    }

    void save_buildings_after_extrusion(const LOD& lod) {
      save_mesh(lod, Intermediate_step::EXTRUDED_BUILDING_BOUNDARIES,
      m_path_bu + "buildings_7_extruded_boundaries");
      save_mesh(lod, Intermediate_step::EXTRUDED_BUILDING_FOOTPRINTS,
      m_path_bu + "buildings_8_extruded_footprints");
    }

    void save_roofs_before_extraction(const LOD& lod) {
      save_points(lod, Intermediate_step::BUILDING_ROOF_POINTS,
      m_path_bu + "buildings_9_roof_points");
      save_mesh(lod, Intermediate_step::APPROXIMATE_BUILDING_BOUNDS,
      m_path_bu + "buildings_10_approximate_bounds");
    }

    void save_roofs_after_extraction(const LOD& lod) {
      save_mesh(lod, Intermediate_step::APPROXIMATE_BUILDING_BOUNDS,
      m_path_bu + "buildings_11_partitioning_input");
      save_mesh(lod, Intermediate_step::BUILDING_PARTITIONING_3,
      m_path_bu + "buildings_12_partitioning_output");
      save_mesh(lod, Intermediate_step::BUILDING_WALLS,
      m_path_bu + "buildings_13_walls");
      save_mesh(lod, Intermediate_step::BUILDING_ROOFS,
      m_path_bu + "buildings_14_roofs");
    }

    // Helpers.
    void save_points(
      const LOD& lod,
      const Intermediate_step step,
      const std::string path) {

      Point_set points;
      Point_inserter<Traits> inserter(points);
      lod.points(
        boost::make_function_output_iterator(inserter),
        step);
      m_saver.export_point_set(points, path);
    }

    void save_polylines(
      const LOD& lod,
      const Intermediate_step step,
      const std::string path) {

      Points_container segments;
      Polyline_inserter<Traits> inserter(segments);
      lod.polylines(
        boost::make_function_output_iterator(inserter),
        step);
      m_saver.export_polylines(segments, path);
    }

    void save_mesh(
      const LOD& lod,
      const Intermediate_step step,
      const std::string path) {

      Points vertices; Indices_container faces; Colors fcolors;
      Polygon_inserter<Traits> inserter(faces, fcolors);
      lod.mesh(
        std::back_inserter(vertices),
        boost::make_function_output_iterator(inserter),
        step);
      m_saver.export_polygon_soup(vertices, faces, fcolors, path);
    }

    void recenter_point_cloud() {

      Point_set point_set(true);
      std::ifstream in(
        "/Users/monet/Documents/lod/results/exp46-meeting/florent/data/bu22/bu22-upsampled-oriented.ply",
        std::ios_base::binary);
      in.precision(20);
      in >> point_set;
      in.close();

      FT x = FT(0), y = FT(0), z = FT(0);
      for (std::size_t i = 0; i < point_set.size(); ++i) {
        const auto& point  = point_set.point(i);
        x += point.x();
        y += point.y();
        z += point.z();
      }

      const FT size = static_cast<FT>(point_set.size());
      x /= size; y /= size; z /= size;

      std::vector<Point_3> points;
      points.reserve(point_set.size());

      std::vector<Vector_3> normals;
      normals.reserve(point_set.size());

      for (std::size_t i = 0; i < point_set.size(); ++i) {
        const auto& point  = point_set.point(i);
        const FT xx = point.x() - x;
        const FT yy = point.y() - y;
        const FT zz = point.z() - z;

        points.push_back(Point_3(xx, yy, zz));
        normals.push_back(point_set.normal(i));
      }

      std::stringstream out;
      out.precision(20);
      out <<
			"ply" 				             +  std::string(_NL_) + ""     <<
			"format ascii 1.0"         +  std::string(_NL_) + ""     <<
      "comment VCGLIB generated" +  std::string(_NL_) + ""     <<
			"element vertex " << points.size() << "" + std::string(_NL_) + "" <<
			"property float x"         +  std::string(_NL_) + ""     <<
			"property float y"         +  std::string(_NL_) + ""     <<
			"property float z"         +  std::string(_NL_) + "" 		 <<
      "property float nx"        +  std::string(_NL_) + ""     <<
			"property float ny"        +  std::string(_NL_) + ""     <<
			"property float nz"        +  std::string(_NL_) + "" 		 <<
      "element face 0"           +  std::string(_NL_) + "" 		 <<
      "property list uchar int vertex_indices" +  std::string(_NL_) + "" <<
			"end_header"               +  std::string(_NL_) + "";

      for (std::size_t i = 0; i < points.size(); ++i)
        out << points[i] << " " << normals[i] << std::endl;

      std::ofstream file(
        "/Users/monet/Documents/lod/results/exp46-meeting/florent/data/bu22/bu22-upsampled-oriented-cheated.ply",
        std::ios_base::out);

      file.precision(20);
      file << out.str();
      file.close();
    }
  }; // Wrapper

} // Levels_of_detail
} // CGAL

#endif // CGAL_LOD_WRAPPER_H
