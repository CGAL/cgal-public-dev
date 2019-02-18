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
#include "Parameters.h"
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
    using Parameters = Parameters<FT>;
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

      m_parameters.update_dependent();

      // Detecting building boundaries.
      m_terminal_parser.add_val_parameter("-alpha_2", m_parameters.alpha_shape_size_2);
      m_terminal_parser.add_val_parameter("-cell_2", m_parameters.grid_cell_width_2);

      m_terminal_parser.add_val_parameter("-rg_search_2", m_parameters.region_growing_search_size_2);
      m_terminal_parser.add_val_parameter("-rg_noise_2", m_parameters.region_growing_noise_level_2);
      m_terminal_parser.add_val_parameter("-rg_angle_2", m_parameters.region_growing_angle_2);
      m_terminal_parser.add_val_parameter("-rg_length_2", m_parameters.region_growing_min_length_2);

      // Computing building footprints.
      m_terminal_parser.add_val_parameter("-kn_width_2", m_parameters.kinetic_min_face_width_2);
      m_terminal_parser.add_val_parameter("-kn_inter_2", m_parameters.kinetic_max_intersections_2);
      m_terminal_parser.add_val_parameter("-bfaces_2", m_parameters.min_faces_per_building_2);

      // Computing tree footprints.
      m_terminal_parser.add_val_parameter("-tr_cell_2", m_parameters.tree_grid_cell_width_2);
      m_terminal_parser.add_val_parameter("-tr_height", m_parameters.min_tree_height);
      m_terminal_parser.add_val_parameter("-tr_radius", m_parameters.min_tree_radius);
      m_terminal_parser.add_val_parameter("-tfaces_2", m_parameters.min_faces_per_tree_2);

      // Extrusion.
      m_terminal_parser.add_val_parameter("-extrusion", m_parameters.extrusion_type);

      // Detecting building roofs. 
      m_terminal_parser.add_val_parameter("-rg_search_3", m_parameters.region_growing_search_size_3);
      m_terminal_parser.add_val_parameter("-rg_noise_3", m_parameters.region_growing_noise_level_3);
      m_terminal_parser.add_val_parameter("-rg_angle_3", m_parameters.region_growing_angle_3);
      m_terminal_parser.add_val_parameter("-rg_area_3", m_parameters.region_growing_min_area_3);

      m_terminal_parser.add_val_parameter("-rc_size", m_parameters.roof_cleaner_min_size);

      // Computing building roofs.
      m_terminal_parser.add_val_parameter("-kn_inter_3", m_parameters.kinetic_max_intersections_3);
      m_terminal_parser.add_val_parameter("-gc_beta_3", m_parameters.graph_cut_beta_3);

      // Info.
      m_parameters.save(m_path);
    }

    void load_input_data() {

      std::cout << std::endl << "Input data: " << std::endl;
      std::ifstream file(m_parameters.data.c_str(), std::ios_base::binary);

      file >> m_point_set;
      file.close();

      std::cout << "File contains " << m_point_set.size() << " points" << std::endl;
        
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

      // Step 1: reconstruct ground as a plane.
      lod.compute_planar_ground();

      Points pg;
      lod.output_ground_as_polygon(std::back_inserter(pg));
      m_saver.export_planar_ground(
        pg, 
        m_path01 + "1_planar_ground");

      // Step 2: detect building boundaries.
      lod.detect_building_boundaries(
        m_parameters.alpha_shape_size_2,
        m_parameters.grid_cell_width_2,
        m_parameters.region_growing_search_size_2,
        m_parameters.region_growing_noise_level_2,
        m_parameters.region_growing_angle_2,
        m_parameters.region_growing_min_length_2);

      Point_set bbpts;
      lod.output_points_along_building_boundary(bbpts.point_back_inserter());
      m_saver.export_point_set(
        bbpts, 
        m_path01 + "2_building_boundary_points");

      Point_set bwpts;
      Insert_point_colored_by_index<Traits> bwp_inserter(bwpts);
      lod.output_points_along_building_walls(
        boost::make_function_output_iterator(bwp_inserter));
      m_saver.export_point_set(
        bwpts, 
        m_path01 + "3_building_wall_points");

      Points_container bbedgs;
      Add_polyline_from_segment<Traits> abbe_adder(bbedgs);
      lod.output_building_boundaries_as_polylines(
        boost::make_function_output_iterator(abbe_adder));
      m_saver.export_polylines(
        bbedgs, 
        m_path01 + "4_approximate_building_boundaries");

      // Step 3: compute building footprints.
      lod.compute_building_footprints(
        m_parameters.kinetic_min_face_width_2,
        m_parameters.kinetic_max_intersections_2,
        m_parameters.min_faces_per_building_2);

      Points vertices; Indices_container faces; Colors fcolors;
      Add_polygon_with_color pr_adder(faces, fcolors);
      
      lod.output_building_partitioning_as_polygon_soup(
        std::back_inserter(vertices),
        boost::make_function_output_iterator(pr_adder));

      m_saver.export_polygon_soup(
        vertices, faces, fcolors,
        m_path01 + "5_building_partitioning");

      vertices.clear(); faces.clear(); fcolors.clear();
      Add_polygon_with_color pr_with_vis_adder(faces, fcolors, true);

      lod.output_building_partitioning_as_polygon_soup(
        std::back_inserter(vertices),
        boost::make_function_output_iterator(pr_with_vis_adder));

      m_saver.export_polygon_soup(
        vertices, faces, fcolors,
        m_path01 + "6_building_visibility");

      vertices.clear(); faces.clear(); fcolors.clear();
      Add_triangle_with_color bfp_adder(faces, fcolors);

      lod.output_building_footprints_as_triangle_soup(
        std::back_inserter(vertices),
        boost::make_function_output_iterator(bfp_adder));

      m_saver.export_polygon_soup(
        vertices, faces, fcolors,
        m_path01 + "7_building_footprints");

      bbedgs.clear();
      Add_polyline_from_segment<Traits> ebbe_adder(bbedgs);
      lod.output_building_boundaries_as_polylines(
        boost::make_function_output_iterator(ebbe_adder));
      m_saver.export_polylines(
        bbedgs, 
        m_path01 + "8_exact_building_boundaries");

      // Step 4: compute tree footprints.
      lod.compute_tree_footprints(
        m_parameters.tree_grid_cell_width_2,
        m_parameters.min_tree_height,
        m_parameters.min_tree_radius,
        m_parameters.min_faces_per_tree_2);

      Point_set trpts;
      Insert_point_colored_by_index<Traits> trp_inserter(trpts);

      lod.output_clustered_vegetation_points(
        boost::make_function_output_iterator(trp_inserter));
      m_saver.export_point_set(
        trpts, 
        m_path01 + "9_tree_points");

      vertices.clear(); faces.clear(); fcolors.clear();
      Add_triangle_with_color tfp_adder(faces, fcolors);

      lod.output_tree_footprints_as_triangle_soup(
        std::back_inserter(vertices),
        boost::make_function_output_iterator(tfp_adder));
      m_saver.export_polygon_soup(
        vertices, faces, fcolors,
        m_path01 + "10_tree_footprints");

      Points_container tredgs;
      Add_polyline_from_segment<Traits> tre_adder(tredgs);
      lod.output_tree_boundaries_as_polylines(
        boost::make_function_output_iterator(tre_adder));
      m_saver.export_polylines(
        tredgs, 
        m_path01 + "11_tree_boundaries");

      // Step 5: LOD0.

      // Step 6: reconstruct smooth ground.
      lod.compute_smooth_ground();

      vertices.clear(); faces.clear(); fcolors.clear();
      Add_triangle_with_color gfp_adder(faces, fcolors);

      lod.output_ground_as_triangle_soup(
        std::back_inserter(vertices),
        boost::make_function_output_iterator(gfp_adder));

      m_saver.export_polygon_soup(
        vertices, faces, fcolors,
        m_path01 + "12_smooth_ground");

      // Step 7: extrude building footprints.
      lod.extrude_building_footprints(
        static_cast<Extrusion_type>(m_parameters.extrusion_type));

      vertices.clear(); faces.clear(); fcolors.clear();
      Add_triangle_with_color ebfp_adder(faces, fcolors);

      lod.output_building_footprints_as_triangle_soup(
        std::back_inserter(vertices),
        boost::make_function_output_iterator(ebfp_adder),
        true);
        
      m_saver.export_polygon_soup(
        vertices, faces, fcolors,
        m_path01 + "13_extruded_building_footprints");

      // Step 8: extrude tree footprints.
      lod.extrude_tree_footprints(
        static_cast<Extrusion_type>(m_parameters.extrusion_type));

      vertices.clear(); faces.clear(); fcolors.clear();
      Add_triangle_with_color etfp_adder(faces, fcolors);

      lod.output_tree_footprints_as_triangle_soup(
        std::back_inserter(vertices),
        boost::make_function_output_iterator(etfp_adder),
        true);

      m_saver.export_polygon_soup(
        vertices, faces, fcolors,
        m_path01 + "14_extruded_tree_footprints");

      // Step 9: LOD1.

      // Step 10: detect building roofs.
      lod.detect_building_roofs(
        m_parameters.region_growing_search_size_3,
        m_parameters.region_growing_noise_level_3,
        m_parameters.region_growing_angle_3,
        m_parameters.region_growing_min_area_3,
        m_parameters.roof_cleaner_min_size);

      Point_set brpts;
      Insert_point_colored_by_index<Traits> brp_inserter(brpts);
      lod.output_points_along_building_roofs(
        boost::make_function_output_iterator(brp_inserter));
      m_saver.export_point_set(
        brpts, 
        m_path2 + "1_building_roof_points");

      vertices.clear(); faces.clear(); fcolors.clear();
      Add_polygon_with_color abr_adder(faces, fcolors);

      lod.output_building_roofs_as_polygon_soup(
        std::back_inserter(vertices),
        boost::make_function_output_iterator(abr_adder));

      m_saver.export_polygon_soup(
        vertices, faces, fcolors,
        m_path2 + "2_approximate_building_roofs");

      // Step 11: compute building roofs.
      lod.compute_building_roofs(
        m_parameters.kinetic_max_intersections_3,
        m_parameters.graph_cut_beta_3);

      vertices.clear(); faces.clear(); fcolors.clear();
      Add_polygon_with_color pi_adder(faces, fcolors);

      lod.output_building_partitioning_in_3_as_polygon_soup(
        std::back_inserter(vertices),
        boost::make_function_output_iterator(pi_adder));

      m_saver.export_polygon_soup(
        vertices, faces, fcolors,
        m_path2 + "3_partitioning_input");

      // Step 12: fit tree icons.

      // Step 13: LOD2.
    }

  }; // Wrapper
    
} // Levels_of_detail
} // CGAL

#endif // CGAL_LOD_WRAPPER_H
