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

      // Fitting tree models.
      m_terminal_parser.add_val_parameter("-tr_precision", m_parameters.tree_precision);

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
    }

  }; // Wrapper
    
} // Levels_of_detail
} // CGAL

#endif // CGAL_LOD_WRAPPER_H
