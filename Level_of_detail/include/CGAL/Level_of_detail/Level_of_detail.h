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
#include <CGAL/Level_of_detail/Reconstruction/Reconstruction_include.h>
#include <CGAL/Level_of_detail/Regularization/Regularization_include.h>
#include <CGAL/Level_of_detail/Shape_detection/Shape_detection_include.h>
#include <CGAL/Level_of_detail/Partitioning/Partitioning_include.h>
#include <CGAL/Level_of_detail/Visibility/Visibility_include.h>
#include <CGAL/Level_of_detail/Buildings/Buildings_include.h>

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
			using Point_3 = typename Kernel::Point_3;
			using Plane_3 = typename Kernel::Plane_3;

			using Parameters 	 = LOD::Parameters<FT>;
			using Data_structure = LOD::Data_structure<Kernel, Input_range, Point_map>;

			using Point_identifier  = typename Data_structure::Point_identifier;
			using Point_identifiers = typename Data_structure::Point_identifiers;

			using Point_map_2 = LOD::Point_from_reference_property_map_2<Point_identifier, Point_2, Point_map>;
			using Point_map_3 = LOD::Dereference_property_map<Point_identifier, Point_map>;

			using Semantic_data_splitter = LOD::Semantic_data_splitter;
			using Visibility_setter 	 = LOD::Visibility_setter;

			using Plane_to_points_fitter = LOD::Plane_to_points_fitter<Kernel>;
			using Bounding_box_estimator = LOD::Bounding_box_estimator<Kernel>;

			using Grid_based_filtering   = LOD::Grid_based_filtering<Kernel, Point_identifier>;
			using Alpha_shapes_filtering = LOD::Alpha_shapes_filtering<Kernel, Point_identifier>;
			
			using Points_tree_2 		       = LOD::Kd_tree_with_data_creator<Kernel, Point_identifier, Point_identifiers, Point_map_2>;
			using Tree_based_lines_estimator_2 = LOD::Tree_based_lines_estimator<Kernel, Point_identifiers, Point_map_2, Points_tree_2>;

			using Linearity_based_sorting_2   = LOD::Scores_based_sorting<Tree_based_lines_estimator_2>;
			using Region_growing_2_normal_map = LOD::Estimated_normal_property_map_2<Kernel, Tree_based_lines_estimator_2>;
			using Region_growing_2       	  = LOD::Points_based_region_growing_2<Kernel, Points_tree_2>;

			using Regularized_segments 		= typename Data_structure::Regularized_segments;
			using Segment_from_region_map_2 = LOD::Segment_from_region_property_map_2<Kernel, Point_identifiers, Point_map_2, Regularized_segments>;
			
			using Segment_regularizer_parameters = LOD::Segment_regularizer_parameters<FT>;
			using Segment_regularizer_2			 = LOD::Segment_regularizer_2<Kernel>;

			using Partition_face_2 			   = typename Data_structure::Partition_face_2;
			using Kinetic_based_partitioning_2 = LOD::Kinetic_based_partitioning_2<Kernel, Partition_face_2>;

			using Partition_point_map 			       = CGAL::Identity_property_map<Point_2>;
			using Triangulation						   = typename Data_structure::Triangulation;
			using Constrained_triangulation_creator    = LOD::Constrained_triangulation_creator<Kernel, Triangulation>;
			using Triangulation_visibility_consistency = LOD::Visibility_consistency<Triangulation>;

			using Building_map 	   = LOD::Building_with_segment_constraints_property_map<Kernel, Triangulation>;
			using Buildings_setter = LOD::Buildings_setter;
			
			using Building 			 = typename Data_structure::Building;
			using Buildings_creator  = LOD::Buildings_creator<Kernel, Building>;
			using Buildings_outliner = LOD::Buildings_outliner<Kernel, Building>;
			
			using Triangulation_building_interior_points_setter = LOD::Buildings_interior_points_setter<Kernel, Triangulation>;
			using Building_height_map 						    = LOD::Building_height_property_map<Kernel, Triangulation, Building>;
			using Buildings_height_setter 		   	   		    = LOD::Buildings_height_setter;

			using Mesh 	= typename Data_structure::Mesh;
			using Lod_0 = typename Data_structure::Lod_0;
			using Lod_1 = typename Data_structure::Lod_1;

			Level_of_detail(const Input_range &input_range, const Point_map &point_map, const Parameters &parameters) :
			m_data_structure(input_range, point_map),
			m_parameters(parameters),
			m_point_map_2(m_data_structure.point_map()),
			m_point_map_3(m_data_structure.point_map()) { 
				
				CGAL_assertion(input_range.size() != 0);
			}

			//////////////////////////////
			// Functions to be documented!

			template<class Semantic_element_map, class Visibility_map_2>
			void build(const Semantic_element_map &semantic_element_map, const Visibility_map_2 &visibility_map_2) {

				if (m_parameters.verbose()) std::cout << std::endl << "... building LOD data ..." << std::endl << std::endl;

				split_semantic_data(semantic_element_map);
				
				fit_ground_plane();
				
				extract_building_boundaries();
				
				simplify_building_boundaries();

				detect_lines();

				regularize_segments();

				create_partitioning();

				compute_visibility(visibility_map_2);

				create_triangulation();

				find_buildings();

				find_building_walls();

				fit_flat_building_roofs();
			}

			template<class Lod>
			void get_lod(Lod &lod) {
				
				if (m_parameters.verbose()) std::cout << "* constructing " << lod.name() << std::endl;
				lod.reconstruct(m_data_structure.buildings(), m_data_structure.ground_bounding_box());
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
				bounding_box_estimator.compute_bounding_box_3(m_data_structure.ground_points(), m_point_map_3, m_data_structure.ground_plane(), m_data_structure.ground_bounding_box());
				
				auto  it = m_data_structure.ground_bounding_box().begin();
				const Point_3 &p1 = *it; ++it;
				const Point_3 &p2 = *it; ++it;
				const Point_3 &p3 = *it; ++it;

				m_data_structure.ground_plane() = Plane_3(p1, p2, p3);
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

				m_data_structure.building_boundary_points().clear();
			}

			void simplify_building_boundaries() {
				if (m_parameters.verbose()) std::cout << "* simplifying building boundary points";

				// Here, we apply grid-based simplification to all building boundary points.
				if (m_parameters.no_simplification()) {
					
					m_data_structure.simplified_building_boundary_points() = m_data_structure.filtered_building_boundary_points();
					if (m_parameters.verbose()) std::cout << " - skipped";

				} else {

					const Grid_based_filtering grid_based_filtering(m_parameters.grid_cell_width());
					grid_based_filtering.apply(m_data_structure.filtered_building_boundary_points(), m_point_map_2, m_data_structure.simplified_building_boundary_points());	
				}

				if (m_parameters.verbose()) std::cout << std::endl;
				m_data_structure.filtered_building_boundary_points().clear();
			}

			void detect_lines() {
				if (m_parameters.verbose()) std::cout << "* detecting lines along building boundaries" << std::endl;

				// In this step, we apply region growing to detect 2D lines that form building boundaries.
				const Points_tree_2 points_tree_2(
					m_data_structure.simplified_building_boundary_points(),
					m_point_map_2, 
					m_parameters.region_growing_2_cluster_epsilon());

				const Tree_based_lines_estimator_2 lines_estimator_2(
					m_data_structure.simplified_building_boundary_points(), 
					m_point_map_2, 
					points_tree_2);
				
				const Linearity_based_sorting_2 linearity_based_sorting_2(lines_estimator_2);
				std::stable_sort(m_data_structure.simplified_building_boundary_points().begin(), m_data_structure.simplified_building_boundary_points().end(), linearity_based_sorting_2);
				
				const Region_growing_2_normal_map normal_map(
					m_data_structure.simplified_building_boundary_points(), 
					m_point_map_2, 
					lines_estimator_2.lines_2());

				Region_growing_2 region_growing_2(
					m_parameters.region_growing_2_epsilon(), 
					m_parameters.region_growing_2_cluster_epsilon(), 
					m_parameters.region_growing_2_normal_threshold(),
					m_parameters.region_growing_2_min_points(),
					points_tree_2);

				region_growing_2.detect(m_data_structure.simplified_building_boundary_points(), m_point_map_2, normal_map, m_data_structure.detected_2d_regions());
				m_data_structure.simplified_building_boundary_points().clear();
			}

			void regularize_segments() {
				if (m_parameters.verbose()) std::cout << "* regularizing segments detected along building boundaries";

				// Here, we regularize segments that form building boundaries wrt to angles and ordinates.
				const Segment_from_region_map_2 segment_from_region_map_2(m_data_structure.detected_2d_regions(), m_point_map_2);

				Segment_regularizer_parameters segment_regularizer_parameters;				
				segment_regularizer_parameters.max_angle_in_degrees() 	  = m_parameters.segment_regularizer_2_max_angle_in_degrees();
				segment_regularizer_parameters.max_difference_in_meters() = m_parameters.segment_regularizer_2_max_difference_in_meters();

				if (m_parameters.no_regularization()) {
					if (m_parameters.verbose()) std::cout << " - skipped ";

					segment_regularizer_parameters.optimize_angles() 	= false;
					segment_regularizer_parameters.optimize_ordinates() = false;
				}
				
				Segment_regularizer_2 segment_regularizer_2(segment_regularizer_parameters);
				segment_regularizer_2.regularize(m_data_structure.detected_2d_regions(), segment_from_region_map_2, m_data_structure.regularized_segments());

				if (m_parameters.verbose()) std::cout << std::endl;
				m_data_structure.detected_2d_regions().clear();
			}

			void create_partitioning() {
				if (m_parameters.verbose()) std::cout << "* computing partitioning" << std::endl;

				// In this step, we create a 2D partitioning of the domain.
				const Kinetic_based_partitioning_2 kinetic_based_partitioning_2(
					m_parameters.kinetic_partitioning_2_num_intersections(),
					m_parameters.kinetic_partitioning_2_min_face_width());
					
				kinetic_based_partitioning_2.compute(m_data_structure.regularized_segments(), m_data_structure.regularized_segment_map(), m_data_structure.partition_faces_2());
			}

			template<class Visibility_map_2>
			void compute_visibility(const Visibility_map_2 &visibility_map_2) {
				if (m_parameters.verbose()) std::cout << "* computing visibility" << std::endl;

				// Here, we try to guess which of the partitioning faces is inside or outside the building.
				const Visibility_setter visibility_setter;
				visibility_setter.set_labels(visibility_map_2, m_data_structure.partition_faces_2());
			}

			void create_triangulation() {
				if (m_parameters.verbose()) std::cout << "* creating triangulation" << std::endl;

				// In this step, we build constrained Delaunay triangulation.
				Partition_point_map partition_point_map;

				const Constrained_triangulation_creator constrained_triangulation_creator;
				constrained_triangulation_creator.make_triangulation_with_info(
					m_data_structure.partition_faces_2(), 
					partition_point_map, 
					m_data_structure.triangulation());

				if (!m_parameters.no_consistent_visibility()) {

					const Triangulation_visibility_consistency triangulation_visibility_consistency;
					triangulation_visibility_consistency.make_consistent(m_data_structure.triangulation());
				}

				m_data_structure.partition_faces_2().clear();
			}

			void find_buildings() {
				if (m_parameters.verbose()) std::cout << "* searching for buildings" << std::endl;

				// Here, we search for sets of triangles that form buildings.
				const Building_map face_to_building_map(
					m_data_structure.triangulation(),
					m_data_structure.regularized_segments(),
					m_data_structure.regularized_segment_map(),
					m_parameters.segment_constraints_threshold());
				
				const Buildings_setter buildings_setter;
				buildings_setter.set_buildings(face_to_building_map, m_data_structure.triangulation());

				const Buildings_creator buildings_creator(m_parameters.min_num_building_floor_faces());
				buildings_creator.create(m_data_structure.triangulation(), m_data_structure.buildings());

				m_data_structure.regularized_segments().clear();
			}

			void find_building_walls() {
				if (m_parameters.verbose()) std::cout << "* searching for building walls" << std::endl;

				// In this step, we search for sets of segments that form building walls.
				const Buildings_outliner buildings_outliner;
				buildings_outliner.find_walls(m_data_structure.triangulation(), m_data_structure.buildings());
			}

			void fit_flat_building_roofs() {
				if (m_parameters.verbose()) std::cout << "* fitting flat building roofs" << std::endl;

				// Here, we fit flat roofs to all buildings with the average (see parameters) building height.
				const Triangulation_building_interior_points_setter triangulation_building_interior_points_setter;
				triangulation_building_interior_points_setter.set_points(m_data_structure.building_interior_points(), m_point_map_2, m_data_structure.triangulation());

				const Building_height_map lod_building_height_map(
					m_data_structure.triangulation(), 
					m_point_map_3,
					m_data_structure.ground_plane(),
					m_parameters.flat_roof_type());

				const Buildings_height_setter buildings_height_setter;
				buildings_height_setter.set_heights(lod_building_height_map, m_data_structure.buildings());

				m_data_structure.building_interior_points().clear();
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