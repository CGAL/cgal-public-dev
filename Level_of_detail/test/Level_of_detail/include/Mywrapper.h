#ifndef CGAL_LEVEL_OF_DETAIL_MYWRAPPER_H
#define CGAL_LEVEL_OF_DETAIL_MYWRAPPER_H

#if defined(WIN32) || defined(_WIN32) 
#define _SR_ "\\"
#else 
#define _SR_ "/" 
#endif

// STL includes.
#include <string>
#include <iostream>

// CGAL includes.
#include <CGAL/Timer.h>
#include <CGAL/property_map.h>

// LOD includes.
#include <CGAL/Level_of_detail_include.h>

// Local includes.
#include "debugging/Mylog.h"
#include "loaders/Myloader.h"
#include "terminal/Myterminal_parser.h"

namespace CGAL {

	namespace Level_of_detail {

		namespace LOD = CGAL::Level_of_detail;

		template<class InputKernel>
		class Mywrapper {

		public:
			using Kernel = InputKernel;
			
			using Parameters = char**;
			using FT 		 = typename Kernel::FT;
			using Point_2 	 = typename Kernel::Point_2;
			using Segment_2  = typename Kernel::Segment_2;
			
			using Log 			  = LOD::Mylog;
			using Loader 		  = LOD::Myloader<Kernel>;
			using Terminal_parser = LOD::Myterminal_parser<FT>;

			using Label = typename Loader::Label;
			
			using LOD_input     = typename Loader::Container;
			using LOD_point_map = typename Loader::Point_map;
			using LOD_label_map = typename Loader::Label_map;

			using LOD_point 	  = typename Loader::Point_3;
			using LOD_point_index = typename Loader::Point_index;

			using LOD_semantic_element_map = LOD::Semantic_element_property_map<LOD_point_index, LOD_label_map>;

			using LOD_base 		 = LOD::Level_of_detail<Kernel, LOD_input, LOD_point_map>;
			using LOD_parameters = typename LOD_base::Parameters;

			using LOD_dereference_point_map = LOD::Dereference_property_map<typename LOD_base::Point_identifier, LOD_point_map>;
			using LOD_partition_point_map   = LOD::Partition_point_property_map<Point_2, LOD_point>;
			
			using LOD_partition_face_2 = LOD::Partition_element<Kernel, CGAL::Polygon_2<Kernel> >;
			using LOD_visibility_map_2 = LOD::Visibility_from_classification_property_map_2<LOD_partition_face_2, Kernel, LOD_input, LOD_point_map, LOD_label_map>;
			
			using LOD_colour_map 			= LOD::Colour_property_map;
			using LOD_visibility_colour_map = LOD::Visibility_colour_property_map;

			Mywrapper(const int num_parameters, const Parameters parameters, const std::string &logs_path) : 
			m_terminal_parser(num_parameters, parameters, logs_path),
			m_logs_path_0_1(logs_path + "tmp" + std::string(_SR_) + "lod_0_1" + std::string(_SR_))
			{ }

			void run_lod_pipeline() {
				
				// First, we set LOD parameters parsing a terminal input or to the default values, including the path to the input data.
				set_lod_parameters();

				// Then, we load data from a file. I use my own custom loader here.
				load_lod_input_data();

				// We also set all necessary property maps.
				set_lod_property_maps();

				// Here, we build all steps of the LOD pipeline with the default template classes.
				// compute_lod_default();
				
				// Here, we build all steps of the LOD pipeline with the user-defined template classes.
				compute_lod_custom();
			}

		private:
			Terminal_parser m_terminal_parser;
			LOD_parameters  m_lod_parameters;
			LOD_input       m_lod_input;

			LOD_point_map m_lod_point_map;
			LOD_label_map m_lod_label_map;

			std::string m_logs_path_0_1;

			void set_lod_parameters() {

				// Set all parameters that can be loaded from the terminal.
				// add_str_parameter  - adds a string-type parameter
				// add_val_parameter  - adds a scalar-type parameter
				// add_bool_parameter - adds a boolean parameter

				std::cout << "Input parameters: " << std::endl;

				m_terminal_parser.add_str_parameter("-data", m_lod_parameters.path_to_input());

				m_terminal_parser.add_bool_parameter("-silent" 	  		 , m_lod_parameters.silent());
				m_terminal_parser.add_bool_parameter("-verbose"	  		 , m_lod_parameters.verbose());
				m_terminal_parser.add_bool_parameter("-no_simplification", m_lod_parameters.no_simplification());
				m_terminal_parser.add_bool_parameter("-no_regularization", m_lod_parameters.no_regularization());
				
				m_terminal_parser.add_val_parameter("-scale", m_lod_parameters.scale());
				m_terminal_parser.add_val_parameter("-eps"  , m_lod_parameters.epsilon());
				
				m_lod_parameters.update_dependent();

				m_terminal_parser.add_val_parameter("-alpha", m_lod_parameters.alpha_shape_size());
				m_terminal_parser.add_val_parameter("-cell" , m_lod_parameters.grid_cell_width());

				m_terminal_parser.add_val_parameter("-rg_nt_2d" , m_lod_parameters.region_growing_2_normal_threshold());
				m_terminal_parser.add_val_parameter("-rg_min_2d", m_lod_parameters.region_growing_2_min_points());
				m_terminal_parser.add_val_parameter("-rg_eps_2d", m_lod_parameters.region_growing_2_epsilon());
				m_terminal_parser.add_val_parameter("-rg_ce_2d" , m_lod_parameters.region_growing_2_cluster_epsilon());

				m_terminal_parser.add_val_parameter("-angle", m_lod_parameters.segment_regularizer_2_max_angle_in_degrees());
			}

			void load_lod_input_data() {

				std::cout << std::endl << "Input data: " << std::endl;
				Loader loader; loader.get_data(m_lod_parameters.path_to_input(), m_lod_input);

				std::cout << "LOD input data are loaded: number of points: " << m_lod_input.number_of_points() << std::endl;
			}

			void set_lod_property_maps() {

				m_lod_point_map = m_lod_input.point_map();
				m_lod_label_map = m_lod_input. template property_map<Label>("label").first;
			}

			void compute_lod_default() {
				
				// Begin.
				CGAL::Timer timer;

				std::cout << std::endl << "Default LOD pipeline... " << std::endl;
				LOD_base lod_base(m_lod_input, m_lod_point_map, m_lod_parameters);
				
				timer.start();

				LOD_semantic_element_map lod_semantic_element_map(m_lod_label_map);
				LOD_visibility_map_2 	 lod_visibility_map_2(m_lod_input, m_lod_point_map, m_lod_label_map);

				lod_base.build(lod_semantic_element_map, lod_visibility_map_2);

				std::cout << std::endl;
				lod_base.get_lod0();
				lod_base.get_lod1();

				// End.
				timer.stop();
				std::cout << std::endl << "Running time: " << timer.time() << " seconds." << std::endl << std::endl;
			}

			void compute_lod_custom() {

				// Begin.
				Log log;
				CGAL::Timer timer;

				std::cout << std::endl << "Custom LOD pipeline... " << std::endl;
				LOD_base lod_base(m_lod_input, m_lod_point_map, m_lod_parameters);

				timer.start();


				// * Step ->
				const LOD_semantic_element_map lod_semantic_element_map(m_lod_label_map);
				lod_base.split_semantic_data(lod_semantic_element_map);
				
				const LOD_dereference_point_map lod_dereference_point_map(m_lod_point_map);
				if (!m_lod_parameters.silent()) {
				
					log.save_points(lod_base.get_internal_data_structure().ground_points(), lod_dereference_point_map, m_logs_path_0_1 + "0_ground_points");
					log.save_points(lod_base.get_internal_data_structure().building_boundary_points(), lod_dereference_point_map, m_logs_path_0_1 + "1_building_boundary_points");
					log.save_points(lod_base.get_internal_data_structure().building_interior_points(), lod_dereference_point_map, m_logs_path_0_1 + "2_building_interior_points");
				}


				// * Step ->
				lod_base.fit_ground_plane();
				if (!m_lod_parameters.silent()) {
					
					Identity_property_map<LOD_point> ground_point_map;
					LOD_colour_map ground_colour_map(LOD::Colour_map_type::RANDOM);

					const std::vector<typename LOD_base::Data_structure::Polygon_3> faces = {{ lod_base.get_internal_data_structure().ground_bounding_box() }};
					log.save_faces(faces, ground_point_map, ground_colour_map, m_logs_path_0_1 + "3_ground_plane");
				}


				// * Step ->
				lod_base.extract_building_boundaries();
				if (!m_lod_parameters.silent()) 
					log.save_points(lod_base.get_internal_data_structure().filtered_building_boundary_points(), lod_dereference_point_map, m_logs_path_0_1 + "4_filtered_building_boundary_points");


				// * Step ->
				lod_base.simplify_building_boundaries();
				if (!m_lod_parameters.silent()) 
					log.save_points(lod_base.get_internal_data_structure().simplified_building_boundary_points(), lod_dereference_point_map, m_logs_path_0_1 + "5_simplified_building_boundary_points");


				// * Step ->
				LOD_colour_map lod_regions_colour_map(LOD::Colour_map_type::RANDOM);

				lod_base.detect_lines();
				if (!m_lod_parameters.silent()) 
					log.save_regions(lod_base.get_internal_data_structure().detected_2d_regions(), lod_dereference_point_map, lod_regions_colour_map, m_logs_path_0_1 + "6_detected_2d_regions");

				
				// * Step ->
				const CGAL::Identity_property_map<Segment_2> segment_map_2;
				
				lod_base.regularize_segments();
				if (!m_lod_parameters.silent())
					log.save_segments(lod_base.get_internal_data_structure().regularized_segments(), segment_map_2, m_logs_path_0_1 + "7_regularized_segments");


				// * Step ->
				const LOD_partition_point_map lod_partition_point_map;
				LOD_colour_map lod_partition_colour_map(LOD::Colour_map_type::RANDOM);

				lod_base.create_partitioning();
				if (!m_lod_parameters.silent())
					log.save_faces(lod_base.get_internal_data_structure().partition_faces_2(), lod_partition_point_map, lod_partition_colour_map, m_logs_path_0_1 + "8_partition");


				// * Step ->
				LOD_visibility_map_2 	  lod_visibility_map_2(m_lod_input, m_lod_point_map, m_lod_label_map);
				LOD_visibility_colour_map lod_visibility_colour_map;

				lod_base.compute_visibility(lod_visibility_map_2);
				if (!m_lod_parameters.silent())
					log.save_faces(lod_base.get_internal_data_structure().partition_faces_2(), lod_partition_point_map, lod_visibility_colour_map, m_logs_path_0_1 + "9_visibility");


				// * Step ->
				std::cout << std::endl;
				lod_base.get_lod0();
				lod_base.get_lod1();


				// End.
				timer.stop();
				std::cout << std::endl << "Running time: " << timer.time() << " seconds." << std::endl << std::endl;
			}
		};
	
	} // Level_of_detail

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_MYWRAPPER_H