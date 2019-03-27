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
    m_path01(m_path + "lod_0_1" + std::string(_SR_)),
    m_path2(m_path + "lod_2" + std::string(_SR_)) 
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
    std::string m_path, m_path01, m_path2;
    Point_set m_point_set;
    Label_map m_label_map;

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

      // Update.
      m_parameters.update_dependent();


      // Clustering buildings.
      m_terminal_parser.add_val_parameter("-bu_clust", m_parameters.buildings.cluster_scale);

      // Detecting building boundaries.
      m_terminal_parser.add_val_parameter("-alpha_2", m_parameters.buildings.alpha_shape_size_2);
      m_terminal_parser.add_val_parameter("-bu_cell_2", m_parameters.buildings.grid_cell_width_2);

      m_terminal_parser.add_val_parameter("-rg_scale_2", m_parameters.buildings.region_growing_scale_2);
      m_terminal_parser.add_val_parameter("-rg_noise_2", m_parameters.buildings.region_growing_noise_level_2);
      m_terminal_parser.add_val_parameter("-rg_angle_2", m_parameters.buildings.region_growing_angle_2);
      m_terminal_parser.add_val_parameter("-rg_length_2", m_parameters.buildings.region_growing_min_length_2);

      // Computing building footprints.
      m_terminal_parser.add_val_parameter("-kn_width_2", m_parameters.buildings.kinetic_min_face_width_2);
      m_terminal_parser.add_val_parameter("-kn_inter_2", m_parameters.buildings.kinetic_max_intersections_2);
      m_terminal_parser.add_val_parameter("-bu_faces_2", m_parameters.buildings.min_faces_per_footprint);
      m_terminal_parser.add_val_parameter("-gc_beta_2", m_parameters.buildings.graph_cut_beta_2);

      // Detecting building roofs. 
      m_terminal_parser.add_val_parameter("-rg_scale_3", m_parameters.buildings.region_growing_scale_3);
      m_terminal_parser.add_val_parameter("-rg_noise_3", m_parameters.buildings.region_growing_noise_level_3);
      m_terminal_parser.add_val_parameter("-rg_angle_3", m_parameters.buildings.region_growing_angle_3);
      m_terminal_parser.add_val_parameter("-rg_area_3", m_parameters.buildings.region_growing_min_area_3);
      m_terminal_parser.add_val_parameter("-roof_scale", m_parameters.buildings.min_roof_scale);

      // Computing building roofs.
      m_terminal_parser.add_val_parameter("-kn_inter_3", m_parameters.buildings.kinetic_max_intersections_3);
      m_terminal_parser.add_val_parameter("-gc_beta_3", m_parameters.buildings.graph_cut_beta_3);

      // Extrusion.
      // m_terminal_parser.add_val_parameter("-bu_extr", m_parameters.buildings.extrusion_type);


      // Clustering trees.
      m_terminal_parser.add_val_parameter("-tr_clust", m_parameters.trees.cluster_scale);

      // Computing tree footprints.
      m_terminal_parser.add_val_parameter("-tr_cell_2", m_parameters.trees.grid_cell_width_2);
      m_terminal_parser.add_val_parameter("-tr_height", m_parameters.trees.min_height);
      m_terminal_parser.add_val_parameter("-tr_radius", m_parameters.trees.min_radius_2);
      m_terminal_parser.add_val_parameter("-tr_faces_2", m_parameters.trees.min_faces_per_footprint);

      // Fitting tree models.
      m_terminal_parser.add_val_parameter("-tr_prec", m_parameters.trees.precision);

      // Extrusion.
      // m_terminal_parser.add_val_parameter("-tr_extr", m_parameters.trees.extrusion_type);


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

      // Define a map from a user-defined label to the LOD semantic label.
      Semantic_map semantic_map(m_label_map, 
      m_parameters.gi,
      m_parameters.bi,
      m_parameters.ii,
      m_parameters.vi);

      // Define a map for computing visibility.
      Visibility_map visibility_map(semantic_map);

      // Create LOD.
      LOD lod(
        m_point_set, 
        m_point_set.point_map(), 
        semantic_map,
        visibility_map);

      std::cout << std::endl << "STEPS:" << std::endl;

      // Ground.
      save_ground(lod, 
      Reconstruction_type::PLANAR_GROUND, m_parameters.ground.precision,
      m_path + "planar_ground");
      save_ground(lod, 
      Reconstruction_type::SMOOTH_GROUND, m_parameters.ground.precision,
      m_path + "smooth_ground");

      // Trees.
      lod.initialize_trees(
        m_parameters.trees.cluster_scale);
      save_tree_clusters(lod);

      lod.compute_tree_footprints(
        m_parameters.trees.grid_cell_width_2,
        m_parameters.trees.min_height,
        m_parameters.trees.min_radius_2,
        m_parameters.trees.min_faces_per_footprint);
      save_trees_before_extrusion(lod);

      lod.extrude_tree_footprints(
        m_parameters.trees.extrusion_type);
      save_trees_after_extrusion(lod);

      // Buildings.

      // LODs.
      save_lod(lod, 
      Reconstruction_type::LOD0, m_parameters.ground.precision, 
      m_path + "LOD0");
      save_lod(lod, 
      Reconstruction_type::LOD1, m_parameters.ground.precision, 
      m_path + "LOD1");
      save_lod(lod, 
      Reconstruction_type::LOD2, m_parameters.ground.precision, 
      m_path + "LOD2");
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

    void save_tree_clusters(const LOD& lod) {
      save_points(lod, Intermediate_step::TREE_CLUSTERS, 
      m_path01 + "trees_1_clusters");
    }

    void save_trees_before_extrusion(const LOD& lod) {
      save_points(lod, Intermediate_step::TREE_POINTS, 
      m_path01 + "trees_2_points");
      save_polylines(lod, Intermediate_step::TREE_BOUNDARIES,
      m_path01 + "trees_3_boundaries");
      save_mesh(lod, Intermediate_step::TREE_FOOTPRINTS,
      m_path01 + "trees_4_footprints");
    }

    void save_trees_after_extrusion(const LOD& lod) {
      save_mesh(lod, Intermediate_step::EXTRUDED_TREE_BOUNDARIES,
      m_path01 + "trees_5_extruded_boundaries");
      save_mesh(lod, Intermediate_step::EXTRUDED_TREE_FOOTPRINTS,
      m_path01 + "trees_5_extruded_footprints");
    }

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
  }; // Wrapper
    
} // Levels_of_detail
} // CGAL

#endif // CGAL_LOD_WRAPPER_H
