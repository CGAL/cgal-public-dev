#ifndef CGAL_LEVEL_OF_DETAIL_H
#define CGAL_LEVEL_OF_DETAIL_H

// STL includes.
#include <iostream>

// LOD includes.
#include <CGAL/Level_of_detail/Parameters.h>
#include <CGAL/Level_of_detail/Enumerations.h>
#include <CGAL/Level_of_detail/Data_structure.h>

#include <CGAL/Level_of_detail/Tools/Data/Semantic_data_splitter.h>
#include <CGAL/Level_of_detail/Tools/Filters/Grid_based_filtering.h>
#include <CGAL/Level_of_detail/Tools/Filters/Alpha_shapes_filtering.h>
#include <CGAL/Level_of_detail/Tools/Fitters/Plane_to_points_fitter.h>
#include <CGAL/Level_of_detail/Tools/General/Bounding_box_estimator.h>

#include <CGAL/Level_of_detail/Tools/Property_maps/Point_property_map_2.h>
#include <CGAL/Level_of_detail/Tools/Property_maps/Dereference_property_map.h>

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

			using Point_identifier      = typename Data_structure::Point_identifier;
			using Dereference_point_map = LOD::Dereference_property_map<Point_identifier, Point_map>;
			using Point_map_2 			= LOD::Point_property_map_2<Point_identifier, Point_2, Point_map>;

			using Semantic_data_splitter = LOD::Semantic_data_splitter<Input_range>;
			using Plane_to_points_fitter = LOD::Plane_to_points_fitter<Kernel>;
			using Bounding_box_estimator = LOD::Bounding_box_estimator<Kernel>;

			using Grid_based_filtering   = LOD::Grid_based_filtering<Kernel, Point_identifier>;
			using Alpha_shapes_filtering = LOD::Alpha_shapes_filtering<Kernel, Point_identifier>;

			Level_of_detail(const Input_range &input_range, const Point_map &point_map, const Parameters &parameters) :
			m_data_structure(input_range, point_map),
			m_parameters(parameters) { 
				
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

				// For the moment we split only ground, building interior, and building boundaries.
				Semantic_data_splitter semantic_data_splitter;
				semantic_data_splitter.split_semantics(m_data_structure.input_range(), semantic_element_map, 
				m_data_structure.ground_points(), m_data_structure.building_boundary_points(), m_data_structure.building_interior_points());
			}

			void fit_ground_plane() {
				if (m_parameters.verbose()) std::cout << "* fitting ground plane" << std::endl;

				// Here, we fit a plane to all ground points.
				CGAL_precondition(m_data_structure.ground_points().size() > 2);
				Dereference_point_map dereference_point_map(m_data_structure.point_map());
				
				Plane_to_points_fitter plane_to_points_fitter;
				plane_to_points_fitter.fit_plane(m_data_structure.ground_points(), dereference_point_map, m_data_structure.ground_plane());

				Bounding_box_estimator bounding_box_estimator;
				bounding_box_estimator.compute_horizontal_bounding_box_3(m_data_structure.ground_points(), dereference_point_map, m_data_structure.ground_plane(), m_data_structure.ground_bounding_box());
			}

			void extract_building_boundaries() {
				if (m_parameters.verbose()) std::cout << "* extracting building boundary points" << std::endl;

				// In this step, we apply alpha shapes to extract only building boundary points.
				Alpha_shapes_filtering alpha_shapes_filtering(m_parameters.alpha());
				m_data_structure.filtered_building_boundary_points().clear();
				
				Point_map_2 point_map_2(m_data_structure.point_map());
				CGAL_precondition(m_data_structure.building_boundary_points().size() > 2 || m_data_structure.building_interior_points().size() > 2);

				if (m_data_structure.building_boundary_points().size() > 2)
					alpha_shapes_filtering.add_points(m_data_structure.building_boundary_points(), point_map_2, m_data_structure.filtered_building_boundary_points());

				if (m_data_structure.building_interior_points().size() > 2)
					alpha_shapes_filtering.add_points(m_data_structure.building_interior_points(), point_map_2, m_data_structure.filtered_building_boundary_points());

				CGAL_postcondition(m_data_structure.filtered_building_boundary_points().size() > 2);
			}

			void simplify_building_boundaries() {
				if (m_parameters.verbose()) std::cout << "* simplifying building boundary points" << std::endl;

				// Here, we apply grid-based simplification to all building boundary points.
				Grid_based_filtering gird_based_filtering;
				m_data_structure.simplified_building_boundary_points().clear();
				
				Point_map_2 point_map_2(m_data_structure.point_map());
				CGAL_precondition(m_data_structure.filtered_building_boundary_points().size() > 0);

				gird_based_filtering.apply(m_data_structure.filtered_building_boundary_points(), point_map_2, m_data_structure.simplified_building_boundary_points());
				CGAL_postcondition(m_data_structure.simplified_building_boundary_points().size() > 0);
			}

			//////////////////////////////////
			// Functions to be not documented!

			inline const Data_structure& get_internal_data_structure() const {
				return m_data_structure;
			}

			// to be added!

		private:
			Data_structure m_data_structure;
			const Parameters &m_parameters;
		};
	
	} // Level_of_detail

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_H