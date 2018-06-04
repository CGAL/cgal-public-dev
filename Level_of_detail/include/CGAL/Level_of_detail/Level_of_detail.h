#ifndef CGAL_LEVEL_OF_DETAIL_H
#define CGAL_LEVEL_OF_DETAIL_H

// STL includes.
#include <memory>
#include <iostream>
#include <algorithm>

// LOD includes.
#include <CGAL/Level_of_detail/Parameters.h>
#include <CGAL/Level_of_detail/Enumerations.h>
#include <CGAL/Level_of_detail/Data_structure.h>

#include <CGAL/Level_of_detail/Tools/Tools_include.h>
#include <CGAL/Level_of_detail/Shape_detection/Shape_detection_include.h>

namespace CGAL {

	namespace Level_of_detail {

		namespace LOD = CGAL::Level_of_detail;

		template<class InputKernel, class InputRange, class InputPointMap>
		class Level_of_detail {

		public:
			using Kernel      = InputKernel;
			using Input_range = InputRange;
			using Point_map   = InputPointMap;

			using FT 	  = typename Kernel::FT;
			using Point_2 = typename Kernel::Point_2;

			using Parameters 	 = LOD::Parameters<FT>;
			using Data_structure = LOD::Data_structure<Kernel, Input_range, Point_map>;

			using Point_identifier  = typename Data_structure::Point_identifier;
			using Point_identifiers = typename Data_structure::Point_identifiers;

			using Point_map_2 = LOD::Point_property_map_2<Point_identifier, Point_2, Point_map>;
			using Point_map_3 = LOD::Dereference_property_map<Point_identifier, Point_map>;

			using Semantic_data_splitter = LOD::Semantic_data_splitter<Input_range>;

			using Plane_to_points_fitter = LOD::Plane_to_points_fitter<Kernel>;
			using Bounding_box_estimator = LOD::Bounding_box_estimator<Kernel>;

			using Grid_based_filtering   = LOD::Grid_based_filtering<Kernel, Point_identifier>;
			using Alpha_shapes_filtering = LOD::Alpha_shapes_filtering<Kernel, Point_identifier>;
			
			using Points_tree_2 		    = LOD::Kd_tree_with_data_creator<Kernel, Point_identifier, Point_identifiers, Point_map_2>;
			using Tree_based_lines_fitter_2 = LOD::Tree_based_lines_fitter<Kernel, Point_identifiers, Point_map_2, Points_tree_2>;

			using Linearity_based_sorting_2   = LOD::Scores_based_sorting<Tree_based_lines_fitter_2>;
			using Region_growing_2_normal_map = LOD::Estimated_normal_property_map_2<Kernel, Tree_based_lines_fitter_2>;
			using Region_growing_2       	  = LOD::Points_based_region_growing_2<Kernel, Points_tree_2>;

			Level_of_detail(const Input_range &input_range, const Point_map &point_map, const Parameters &parameters) :
			m_data_structure(input_range, point_map),
			m_parameters(parameters),
			m_point_map_2(m_data_structure.point_map()),
			m_point_map_3(m_data_structure.point_map()) { 
				
				CGAL_assertion(input_range.size() != 0);
			}

			//////////////////////////////
			// Functions to be documented!

			template<class Semantic_element_map>
			void build(const Semantic_element_map &semantic_element_map) {

				if (m_parameters.verbose()) std::cout << std::endl << "... building LOD data ..." << std::endl << std::endl;

				split_semantic_data(semantic_element_map);
				
				fit_ground_plane();
				
				extract_building_boundaries();
				
				simplify_building_boundaries();

				detect_lines();
			}

			void get_lod0() {
				if (m_parameters.verbose()) std::cout << "* constructing LOD0" << std::endl;
			}

			void get_lod1() {
				if (m_parameters.verbose()) std::cout << "* constructing LOD1" << std::endl;
			}

			template<class Semantic_element_map>
			void split_semantic_data(const Semantic_element_map &semantic_element_map) {
				if (m_parameters.verbose()) std::cout << "* splitting semantic data" << std::endl;

				// In this step, we split only ground, building interior, and building boundaries.
				const Semantic_data_splitter semantic_data_splitter;
				semantic_data_splitter.split_semantics(m_data_structure.input_range(), semantic_element_map, 
				m_data_structure.ground_points(), m_data_structure.building_boundary_points(), m_data_structure.building_interior_points());
			}

			void fit_ground_plane() {
				if (m_parameters.verbose()) std::cout << "* fitting ground plane" << std::endl;

				// Here, we fit a plane to all ground points.
				const Plane_to_points_fitter plane_to_points_fitter;
				plane_to_points_fitter.fit_plane(m_data_structure.ground_points(), m_point_map_3, m_data_structure.ground_plane());

				const Bounding_box_estimator bounding_box_estimator;
				bounding_box_estimator.compute_horizontal_bounding_box_3(m_data_structure.ground_points(), m_point_map_3, m_data_structure.ground_plane(), m_data_structure.ground_bounding_box());

				m_data_structure.ground_points().clear();
			}

			void extract_building_boundaries() {
				if (m_parameters.verbose()) std::cout << "* extracting building boundary points" << std::endl;

				// In this step, we apply alpha shapes to extract only building boundary points.
				CGAL_precondition(m_data_structure.building_boundary_points().size() > 2 || m_data_structure.building_interior_points().size() > 2);
				
				m_data_structure.filtered_building_boundary_points().clear();
				const Alpha_shapes_filtering alpha_shapes_filtering(m_parameters.alpha_shape_size());
				
				if (m_data_structure.building_boundary_points().size() > 2)
					alpha_shapes_filtering.add_points(m_data_structure.building_boundary_points(), m_point_map_2, m_data_structure.filtered_building_boundary_points());

				if (m_data_structure.building_interior_points().size() > 2)
					alpha_shapes_filtering.add_points(m_data_structure.building_interior_points(), m_point_map_2, m_data_structure.filtered_building_boundary_points());
			}

			void simplify_building_boundaries() {
				if (m_parameters.verbose()) std::cout << "* simplifying building boundary points" << std::endl;

				// Here, we apply grid-based simplification to all building boundary points.
				const Grid_based_filtering grid_based_filtering(m_parameters.grid_cell_width());
				grid_based_filtering.apply(m_data_structure.filtered_building_boundary_points(), m_point_map_2, m_data_structure.simplified_building_boundary_points());

				m_data_structure.filtered_building_boundary_points().clear();
			}

			void detect_lines() {
				if (m_parameters.verbose()) std::cout << "* detecting lines along building boundaries" << std::endl;

				// In this step, we apply region growing to detect 2D lines that form building boundaries.
				const Points_tree_2 points_tree_2(
					m_data_structure.simplified_building_boundary_points(),
					m_point_map_2, 
					m_parameters.region_growing_2_cluster_epsilon());

				const Tree_based_lines_fitter_2 lines_fitter_2(
					m_data_structure.simplified_building_boundary_points(), 
					m_point_map_2, 
					points_tree_2);
				
				const Linearity_based_sorting_2 linearity_based_sorting_2(lines_fitter_2);
				std::stable_sort(m_data_structure.simplified_building_boundary_points().begin(), m_data_structure.simplified_building_boundary_points().end(), linearity_based_sorting_2);
				
				const Region_growing_2_normal_map normal_map(
					m_data_structure.simplified_building_boundary_points(), 
					m_point_map_2, 
					lines_fitter_2.lines_2());

				Region_growing_2 region_growing_2(
					m_parameters.region_growing_2_epsilon(), 
					m_parameters.region_growing_2_cluster_epsilon(), 
					m_parameters.region_growing_2_normal_threshold(),
					m_parameters.region_growing_2_min_points(),
					points_tree_2);

				region_growing_2.detect(m_data_structure.simplified_building_boundary_points(), m_point_map_2, normal_map, m_data_structure.detected_2d_regions());
				m_data_structure.simplified_building_boundary_points().clear();
			}

			//////////////////////////////////
			// Functions to be not documented!

			inline const Data_structure& get_internal_data_structure() const {
				return m_data_structure;
			}

		private:
			Data_structure m_data_structure;
			const Parameters &m_parameters;

			const Point_map_2 m_point_map_2;
			const Point_map_3 m_point_map_3;
		};
	
	} // Level_of_detail

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_H