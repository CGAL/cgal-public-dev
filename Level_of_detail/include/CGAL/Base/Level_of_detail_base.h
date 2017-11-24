#ifndef CGAL_LEVEL_OF_DETAIL_BASE_H
#define CGAL_LEVEL_OF_DETAIL_BASE_H

// STL includes.
#include <map>
#include <memory>
#include <string>
#include <iostream>
#include <vector>

// New CGAL includes.
#include <CGAL/Level_of_detail_enum.h>
#include <CGAL/Mylog/Mylog.h>

namespace CGAL {

	namespace LOD {

		// Facade class that accumulates all necessary objects and operations
		// related to the level of detail (LOD) reconstruction.
		template<class LodTraits>
		class Level_of_detail_base {

		public:
			// Main typedefs.
			typedef LodTraits 				      Traits;
			typedef typename Traits::Kernel       Kernel;

			typedef typename Traits::Container_2D Container_2D;
			typedef typename Traits::Container_3D Container_3D;

			typedef typename Traits::Loader       Loader;
			typedef typename Traits::Preprocessor Preprocessor;

			typedef typename Kernel::FT 	   FT;
			typedef typename Kernel::Point_2   Point_2;
			typedef typename Kernel::Point_3   Point_3;
			typedef typename Kernel::Plane_3   Plane_3;
			typedef typename Kernel::Line_2    Line_2;
			typedef typename Kernel::Segment_2 Segment_2;

			typedef typename Traits::Building_boundary_selector Building_boundary_selector; // Maybe use a factory here? 
			typedef typename Traits::Building_interior_selector Building_interior_selector;
			typedef typename Traits::Clutter_selector 		    Clutter_selector;
			typedef typename Traits::Ground_selector 		    Ground_selector;

			typedef typename Traits::Vertical_regularizer Vertical_regularizer;
			typedef typename Traits::Ground_projector     Ground_projector;
			typedef typename Traits::Projected_points     Projected_points;
			typedef typename Traits::Planes        		  Planes;

			typedef typename Traits::Building_splitter 	  Building_splitter;
			typedef typename Traits::Building_outliner 	  Building_outliner;

			typedef typename Traits::Building_min_roof_fitter Building_min_roof_fitter;
			typedef typename Traits::Building_avg_roof_fitter Building_avg_roof_fitter;
			typedef typename Traits::Building_max_roof_fitter Building_max_roof_fitter;

			typedef Planes Boundary_data;

			typedef typename Traits::Structuring_2 	  Structuring_2;
			typedef typename Traits::Visibility_2  	  Visibility_2;
			typedef typename Traits::Region_growing_2 Region_growing_2;
			
			typedef typename Traits::Utils Utils;
			
			typedef typename Traits::CDT        		CDT;
			typedef typename CDT::Vertex_handle 		Vertex_handle;
			typedef typename CDT::Face_handle   		Face_handle;
			typedef typename CDT::Finite_edges_iterator Edge_iterator;
			typedef typename CDT::Finite_faces_iterator Face_iterator;

			typedef typename Traits::Graph_cut Graph_cut;
			typedef typename Traits::Lods Lods;

			typedef typename Traits::Mesh 			   Mesh;
			typedef typename Traits::Mesh_facet_colors Mesh_facet_colors;

			typedef typename Traits::Buildings Buildings;

			typedef typename Lods::Point  Ground_point;
			typedef typename Lods::Ground Ground;

			typedef typename Traits::Grid_simplifier Grid_simplifier;
			typedef typename Traits::Thinning 		 Thinning;

			typedef typename Traits::Clutter_processor Clutter_processor;
			typedef Thinning_fitter_type    		   Clutter_fitter_type;
			typedef Grid_new_point_type 			   Clutter_new_point_type;


			// Extra.
			using Plane_iterator = typename Planes::const_iterator;

			using Index   = int;
			using Indices = std::vector<Index>;

			using Structured_points  = std::vector< std::vector<Point_2> >; 			  
			using Structured_labels  = std::vector< std::vector<Structured_label> >;  
			using Structured_anchors = std::vector< std::vector<std::vector<int> > >;
			
			using Log = CGAL::LOD::Mylog;

			using Lines    = std::vector<Line_2>;
			using Segments = std::vector<Segment_2>;

			using Label     = typename Traits::Label;
			using Label_map = typename Container_3D:: template Property_map<Label>;

			using Point_index = typename Container_3D::Index;
			using Face_points_map = std::map<Face_handle, std::vector<Point_index> >;

			enum class Program_version  { VER0 };
			enum class Pipeline_version { WITH_SHAPE_DETECTION, WITHOUT_SHAPE_DETECTION };


			//////////////
			// Main class!
			Level_of_detail_base() :
			m_prefix_path("/Users/danisimo/Documents/pipeline/data/"),
			m_default_path("default"),
			m_preprocessor_scale(-FT(1)) ,
			m_structuring_epsilon(-FT(1)),
			m_structuring_log(false),
			m_structuring_resample(true),
			m_structuring_get_all_points(false),
			m_add_cdt_clutter(true),
			m_add_cdt_bbox(false),
			m_visibility_save_info(false),
			m_visibility_approach(Visibility_approach::POINT_BASED),
			m_visibility_method(Visibility_method::POINT_BASED_CLASSIFICATION),
			m_visibility_num_samples(0),
			m_graph_cut_save_info(false),
			m_graph_cut_alpha(-FT(1)),
			m_graph_cut_beta(-FT(1)),
			m_graph_cut_gamma(-FT(1)),
			m_building_boundaries_save_internal_info(false), 
			m_building_boundaries_max_inner_iters(0),
			m_building_boundaries_max_outer_iters(0),
			m_roof_fitter_type(Roof_fitter_type::MAX), 
			m_clean_projected_points(true),
			m_max_reg_angle(-FT(1)),
			m_regularizer_reject_planes(true), 
			m_use_boundaries(true),
			m_prog_version(Program_version::VER0),
			m_pipeline_version(Pipeline_version::WITHOUT_SHAPE_DETECTION),
			m_visibility_show_progress(false),
			m_visibility_norm_threshold(-FT(1)),
			m_clutter_knn(0),
			m_clutter_cell_length(-FT(1)),
			m_clutter_fitter_type(Clutter_fitter_type::LINE),
			m_clutter_new_point_type(Clutter_new_point_type::BARYCENTRE),
			m_visibility_num_neighbours(0),
			m_visibility_sampler(Visibility_sampler::RANDOM_UNIFORM_0),
			m_visibility_rays_per_side(0),
			m_visibility_small_edge_threshold(FT(0)),
			m_building_boundary_type(Building_boundary_type::UNORIENTED),
			m_visibility_angle_eps(-FT(1)),
			m_thinning_neighbour_search_type(Neighbour_search_type::KNN),
			m_thinning_fuzzy_radius(-FT(1)),
			m_thinning_type(Thinning_type::NAIVE),
			m_region_growing_epsilon(-FT(1)),
			m_region_growing_cluster_epsilon(-FT(1)),
			m_region_growing_normal_threshold(-FT(1)),
			m_region_growing_min_points(0),
			m_with_region_growing(true),
			m_use_grid_simplifier_first(false),
			m_alpha_shape_size(-FT(1)),
			m_use_alpha_shapes(false),
			m_structuring_corner_algorithm(Structuring_corner_algorithm::GRAPH_BASED),
			m_structuring_adjacency_method(Structuring_adjacency_threshold_method::LOCAL),
			m_structuring_adjacency_value(-FT(1)),
			m_structuring_global_everywhere(true)
			{ } // Do I need to create an instance of these traits here?


			//////////////////
			// Parameter functions!
			void set_prefix_path(const std::string &path) {

				assert(path != "path_to_the_data_folder");
				m_prefix_path = path;
			}


			// Clutter.
			void add_clutter(const bool new_state) {
				m_add_cdt_clutter = new_state;
			}

			void set_clutter_cell_side_length(const FT new_value) {

				assert(new_value > FT(0));
				m_clutter_cell_length = new_value;
			}


			// Region growing.
			void set_region_growing_epsilon(const FT new_value) {

				assert(new_value > FT(0));
				m_region_growing_epsilon = new_value;
			}

			void set_region_growing_cluster_epsilon(const FT new_value) {

				assert(new_value > FT(0));
				m_region_growing_cluster_epsilon = new_value;
			}

			void set_region_growing_normal_threshold(const FT new_value) {

				assert(new_value > FT(0) && new_value < FT(1));
				m_region_growing_normal_threshold = new_value;
			}

			void set_region_growing_min_points(const size_t new_value) {

				assert(new_value > 1);
				m_region_growing_min_points = new_value;
			}


			// Structuring.
			void set_structuring_epsilon(const FT new_value) {

				assert(new_value > FT(0));
				m_structuring_epsilon = new_value;
			}

			void set_structuring_adjacency_value(const FT new_value) {

				assert(new_value > FT(0));
				m_structuring_adjacency_value = new_value;
			}

			void get_all_structuring_points(const bool new_state) {
				m_structuring_get_all_points = new_state;
			}


			// Graph cut.
			void set_graph_cut_beta(const FT new_value) {

				assert(new_value >= FT(0));
				m_graph_cut_beta = new_value;
			}

			void set_graph_cut_gamma(const FT new_value) {

				assert(new_value >= FT(0));
				m_graph_cut_gamma = new_value;
			}


			//////////////////
			// Main functions!

			// Set default parameters.
			void set_default_parameters() {
				
				set_global_parameters();
				assert_global_parameters();
			}


			// All versions.
			void create_lods() {

				switch (m_prog_version) {

					case Program_version::VER0:
						create_lods_ver0();
						break;

					default:
						assert(!"Wrong program version!");
						break;
				}
			}

		private:
			void start_execution(
				Log &log) {

				// Create log and set default parameters.
				std::cout << "\nstarting ..." << std::endl;
				log.out << "START EXECUTION\n\n\n";
			}

			void loading_data(
				Container_3D &input, 
				Log &log, 
				const size_t exec_step) {

				// Read data.
				std::cout << "(" << exec_step << ") loading" << std::endl;

				m_loader.get_data(m_default_path + ".ply", input);

				log.out << "(" << exec_step << ") Data are loaded. Number of points: " << input.number_of_points() << std::endl << std::endl;

				// Log mock_saver; mock_saver.save_ply<Traits, Container_3D>(input, "basic_mock", true);
			}

			void getting_all_planes(
				Planes &all_planes, 
				Log &log, 
				const Container_3D &input, const size_t exec_step) {

				// Find a set of planes related to the points. Basically here we emulate RANSAC.
				// For each plane we store indices of all points contained in this plane.
				std::cout << "(" << exec_step << ") planes" << std::endl;

				const auto number_of_planes = m_preprocessor.get_planes(input, all_planes);
				assert(number_of_planes >= 0);

				log.out << "(" << exec_step << ") Planes are found. Number of planes: " << number_of_planes << std::endl << std::endl;
			}

			void applying_selection(
				Indices &ground_idxs,
				Indices &building_boundary_idxs,
				Indices &building_interior_idxs, 
				Log &log, 
				const Container_3D &input, const size_t exec_step) {

				// Split data with respect to 2 different semantic labels.
				std::cout << "(" << exec_step << ") selection" << std::endl;

				m_ground_selector.select_elements(input, std::back_inserter(ground_idxs));
				m_building_boundary_selector.select_elements(input, std::back_inserter(building_boundary_idxs));
				m_building_interior_selector.select_elements(input, std::back_inserter(building_interior_idxs));

				log.out << "(" << exec_step << " a) Ground is found. Number of elements: " 			   << ground_idxs.size() << std::endl;
				log.out << "(" << exec_step << " b) Building boundaries are found. Number of elements: " << building_boundary_idxs.size() << std::endl << std::endl;
			}

			void ground_fitting(
				Plane_3 &base_ground_plane, 
				Plane_3 &fitted_ground_plane, 
				Log &log, 
				const Indices &ground_idxs, const Container_3D &input, const size_t exec_step) {

				// Create plane from the ground points.
				std::cout << "(" << exec_step << ") ground plane fitting" << std::endl;

				base_ground_plane = Plane_3(FT(0), FT(0), FT(1), FT(0)); // use XY plane instead
				m_utils.fit_ground_plane(input, ground_idxs, fitted_ground_plane);
				
				log.out << "(" << exec_step << " a) Base ground plane is: "        << base_ground_plane   << std::endl;
				log.out << "(" << exec_step << " b) Data-fitted ground plane is: " << fitted_ground_plane << std::endl << std::endl;
			}

			void getting_boundary_points(
				Boundary_data &building_boundaries, 
				Boundary_data &boundary_clutter, 
				Log &log, 
				const Indices &building_boundary_idxs, const Indices &building_interior_idxs, const Container_3D &input, const size_t exec_step) {

				// Map indices from all detected planes to the ones that are a part of the given facades.
				std::cout << "(" << exec_step << ") getting boundaries" << std::endl;
				
				bool with_shape_detection = false;
				if (m_pipeline_version == Pipeline_version::WITH_SHAPE_DETECTION) with_shape_detection = true;

				m_preprocessor.use_alpha_shapes(m_use_alpha_shapes);
				m_preprocessor.set_alpha(m_alpha_shape_size);

				const auto number_of_boundaries = 
				m_preprocessor.get_boundary_points(input, building_boundary_idxs, building_interior_idxs, with_shape_detection, building_boundaries, boundary_clutter);

				log.out << "(" << exec_step << " a) Planes for building's boundary are found. Number of planes: " << number_of_boundaries 		   << std::endl;
				log.out << "(" << exec_step << " b) Boundary clutter is found. Number of points: " 				  << boundary_clutter.at(0).size() << std::endl << std::endl;
			}

			void regularizing(
				Boundary_data &building_boundaries, 
				Container_3D &input,
				Log &log, 
				const Plane_3 &base_ground_plane, const size_t exec_step) {

				// Make all nearly vertical planes in the building's boundary exactly vertical.
				std::cout << "(" << exec_step << ") regularizing" << std::endl;

				m_vertical_regularizer.set_max_regularization_angle(m_max_reg_angle);
				m_vertical_regularizer.reject_planes(m_regularizer_reject_planes);

				const auto number_of_boundaries = building_boundaries.size();
				const auto number_of_regularized_planes = m_vertical_regularizer.regularize(building_boundaries, input, base_ground_plane);

				log.out << "(" << exec_step << ") Building's nearly vertical planes are regularized. Number of regularized planes: " << number_of_regularized_planes <<
				", number of rejected planes: " << number_of_boundaries - building_boundaries.size() << std::endl << std::endl;

				// Log ply_saver; ply_saver.save_ply<Kernel, Container_3D>(input, "regularized", true);
			}

			void projecting(
				Projected_points &building_boundaries_projected, 
				Projected_points &boundary_clutter_projected,
				Log &log, 
				const Plane_3 &base_ground_plane, const Boundary_data &building_boundaries, const Boundary_data &boundary_clutter, const Container_3D &input, const size_t exec_step) {


				// Project all vertical building's boundaries onto the ground plane.
				std::cout << "(" << exec_step << ") projecting; ";

				auto number_of_projected_points = m_ground_projector.project(input, building_boundaries, base_ground_plane, building_boundaries_projected);
				log.out << "(" << exec_step << " a) Building's boundary planar points are projected. Number of projected points: " << number_of_projected_points << std::endl;

				std::cout << "boundaries projected: " << number_of_projected_points << "; ";

				Log points_exporter; 
				if (!building_boundaries_projected.empty()) points_exporter.export_projected_points_as_xyz("tmp/projected_boundaries", building_boundaries_projected, m_default_path);


				// Clutter.
				if (m_pipeline_version == Pipeline_version::WITHOUT_SHAPE_DETECTION || 
				   (m_pipeline_version == Pipeline_version::WITH_SHAPE_DETECTION && m_add_cdt_clutter)) {
					
					number_of_projected_points = m_ground_projector.project(input, boundary_clutter, base_ground_plane, boundary_clutter_projected);
					log.out << "(" << exec_step << " b) Building's boundary clutter is projected. Number of projected points: " << number_of_projected_points << std::endl;
				
					std::cout << "clutter projected: " << number_of_projected_points << "; ";

					points_exporter.clear(); 
					if (!boundary_clutter_projected.empty()) points_exporter.export_projected_points_as_xyz("tmp/projected_clutter", boundary_clutter_projected, m_default_path);
				}

				log.out << std::endl;
				std::cout << std::endl;
			}

			void applying_thinning(
				Projected_points &boundary_clutter_projected, 
				Log &log, 
				const Container_3D &input,
				const size_t exec_step) {

				// Regularize/thin points in the clutter.
				std::cout << "(" << exec_step << ") applying thinning; ";

				m_thinnning.set_number_of_neighbours(m_clutter_knn);
				m_thinnning.set_fitter_type(m_clutter_fitter_type);
				m_thinnning.set_neighbour_search_type(m_thinning_neighbour_search_type);
				m_thinnning.set_fuzzy_radius(m_thinning_fuzzy_radius);
				m_thinnning.set_thinning_type(m_thinning_type);

				Boundary_data stub;
				const auto number_of_thinned_points = m_thinnning.process(stub, boundary_clutter_projected, input);

				std::cout << "thinned points: " << number_of_thinned_points << std::endl;
				log.out << "(" << exec_step << ") Projected points are thinned. Number of thinned points: " << number_of_thinned_points << std::endl << std::endl;
			}

			void applying_grid_simplification(
				Projected_points &boundary_clutter_projected,
				Log &log,  
				const size_t exec_step) {

				// Remove unnecessary points from the clutter.
				std::cout << "(" << exec_step << ") applying grid simplification; ";

				m_grid_simplifier.set_grid_cell_length(m_clutter_cell_length);
				m_grid_simplifier.set_new_point_type(m_clutter_new_point_type);

				Boundary_data stub;
				const auto number_of_removed_points = m_grid_simplifier.process(stub, boundary_clutter_projected);

				std::cout << "removed points: " << number_of_removed_points << std::endl;
				log.out << "(" << exec_step << ") Projected points are simplified. Number of removed points: " << number_of_removed_points << std::endl << std::endl;
			}

			void detecting_2d_lines(
				Boundary_data &boundary_clutter   , Projected_points &boundary_clutter_projected, 
				Boundary_data &building_boundaries, Projected_points &building_boundaries_projected, 
				Log &log, 
				const Container_3D &input,
				const size_t exec_step) {

				// Detect lines in 2D using region growing.
				std::cout << "(" << exec_step << ") detecting 2d lines; ";

				m_region_growing.set_epsilon(m_region_growing_epsilon);
				m_region_growing.set_cluster_epsilon(m_region_growing_cluster_epsilon);
				m_region_growing.set_normal_threshold(m_region_growing_normal_threshold);
				m_region_growing.set_minimum_shape_points(m_region_growing_min_points);

				const auto number_of_detected_lines = m_region_growing.detect(
					boundary_clutter   , boundary_clutter_projected,
					building_boundaries, building_boundaries_projected,
					input);

				std::cout << "detected lines: " << number_of_detected_lines << std::endl;
				log.out << "(" << exec_step << ") 2D lines are detected. Number of detected lines: " << number_of_detected_lines << std::endl << std::endl;
			}

			void cleaning_projected_points(
				Projected_points &building_boundaries_projected, 
				Boundary_data &building_boundaries, 
				Log &log, 
				const size_t exec_step) {

				assert(m_clean_projected_points);

				// Clean projected points by removing all points that lie far away from the center cluster of points.
				std::cout << "(" << exec_step << ") cleaning" << std::endl;
				m_preprocessor.set_scale(m_preprocessor_scale);
					
				const auto number_of_removed_points = m_preprocessor.clean_projected_points(building_boundaries_projected, building_boundaries);
				log.out << "(" << exec_step << ") Building's boundaries are cleaned. Number of removed points: " << number_of_removed_points << std::endl << std::endl;
			}

			void line_fitting(
				Lines &lines, 
				Log &log, 
				const Boundary_data &building_boundaries, const Projected_points &building_boundaries_projected, const size_t exec_step) {

				// Fit lines to the projected points in 2D.
				std::cout << "(" << exec_step << ") line fitting" << std::endl;

				const auto number_of_fitted_lines = m_utils.fit_lines_to_projected_points(building_boundaries_projected, building_boundaries, lines);

				log.out << "(" <<  exec_step << ") Lines are fitted. Number of fitted lines: " << number_of_fitted_lines << std::endl << std::endl;
			}

			void creating_segments(
				Segments &segments, 
				Log &log, 
				const Lines &lines, const Boundary_data &building_boundaries, const Projected_points &building_boundaries_projected, const size_t exec_step) {

				// Find segments from the given lines.
				std::cout << "(" << exec_step << ") creating segments" << std::endl;

				const auto number_of_segments = m_utils.create_segments_from_lines(building_boundaries_projected, building_boundaries, lines, segments);

				log.out << "(" << exec_step << ") Segments are created. Number of created segments: " << number_of_segments << std::endl << std::endl;
				Log segments_exporter; segments_exporter.export_segments_as_obj("tmp/segments", segments, m_default_path);
			}

			void applying_2d_structuring(
				Log &log, 
				const Lines &lines, const Boundary_data &building_boundaries, const Projected_points &building_boundaries_projected, const size_t exec_step) {

				// Apply 2D structuring algorithm.
				std::cout << "(" << exec_step << ") 2d structuring" << std::endl;

				m_structuring = std::make_unique<Structuring_2>(building_boundaries_projected, building_boundaries, lines);
				
				m_structuring->set_epsilon(m_structuring_epsilon);
				m_structuring->save_log(m_structuring_log);
				m_structuring->resample(m_structuring_resample);
				m_structuring->set_corner_algorithm(m_structuring_corner_algorithm);
				m_structuring->set_adjacency_threshold_method(m_structuring_adjacency_method);
				m_structuring->set_adjacency_threshold(m_structuring_adjacency_value);
				m_structuring->set_global_everywhere(m_structuring_global_everywhere);

				const auto number_of_structured_segments = m_structuring->structure_point_set();

				log.out << "(" << exec_step << ") 2D Structuring is applied. Number of structured segments: " << number_of_structured_segments << std::endl << std::endl;
			}

			void processing_clutter(
				Boundary_data &boundary_clutter, 
				Projected_points &boundary_clutter_projected, 
				Log &log, 
				const Container_3D &input,
				const size_t exec_step) {

				// Regularize and remove unnecessary points from the clutter.
				std::cout << "(" << exec_step << ") processing clutter; ";

				m_clutter_processor.set_number_of_neighbours(m_clutter_knn);
				m_clutter_processor.set_grid_cell_length(m_clutter_cell_length);
				m_clutter_processor.set_fitter_type(m_clutter_fitter_type);
				m_clutter_processor.set_new_point_type(m_clutter_new_point_type);
				m_clutter_processor.set_neighbour_search_type(m_thinning_neighbour_search_type);
				m_clutter_processor.set_fuzzy_radius(m_thinning_fuzzy_radius);
				m_clutter_processor.set_thinning_type(m_thinning_type);

				const auto number_of_removed_points = m_clutter_processor.process(boundary_clutter, boundary_clutter_projected, input);

				std::cout << "removed points: " << number_of_removed_points << std::endl;
				log.out << "(" << exec_step << ") Clutter is processed. Number of removed points: " << number_of_removed_points << std::endl << std::endl;
			}

			void creating_cdt(
				CDT &cdt, 
				Log &log, 
				const Boundary_data &boundary_clutter, const Projected_points &boundary_clutter_projected, const Container_3D &input, const size_t exec_step) {

				// Compute constrained Delaunay triangulation of the structured points.
				std::cout << "(" << exec_step << ") creating cdt" << std::endl;

				auto number_of_faces = -1;
				if (m_structuring != nullptr && !m_structuring->is_empty()) {
					
					const Structured_points   &structured_points = m_structuring_get_all_points ?  m_structuring->get_structured_points() :  m_structuring->get_segment_end_points();
					const Structured_labels   &structured_labels = m_structuring_get_all_points ?  m_structuring->get_structured_labels() :  m_structuring->get_segment_end_labels();
					const Structured_anchors &structured_anchors = m_structuring_get_all_points ? m_structuring->get_structured_anchors() : m_structuring->get_segment_end_anchors();

					number_of_faces = m_utils.compute_cdt(structured_points, structured_labels, structured_anchors, m_structuring_adjacency_value, cdt, 
													 m_add_cdt_clutter, boundary_clutter, boundary_clutter_projected, 
													 m_add_cdt_bbox, input);
				} else {

					number_of_faces = m_utils.compute_cdt(Structured_points(), Structured_labels(), Structured_anchors(), m_structuring_adjacency_value, cdt, 
													 m_add_cdt_clutter, boundary_clutter, boundary_clutter_projected, 
													 m_add_cdt_bbox, input);
				}

				assert(number_of_faces != -1);
				log.out << "(" << exec_step << ") Constrained Delaunay triangulation of the structured points is built. Number of faces: " << number_of_faces << std::endl << std::endl;
				assert(!m_add_cdt_bbox); // visibility and graph cut do not work if bbox vertices are added to CDT!
			}

			void converting_3d_to_2d(
				Container_2D &input_2d, 
				Face_points_map &fp_map,
				Log &log, 
				const CDT &cdt, const Container_3D &input, const size_t exec_step) {

				// Convert 3D input to 2D input.			
				std::cout << "(" << exec_step << ") converting 3d input into 2d input and setting face to points map" << std::endl;

				const auto number_of_converted_points = m_utils.get_2d_input_and_face_points_map(cdt, input, input_2d, fp_map);

				log.out << "(" << exec_step << ") 3D input is converted into 2D input and face to points map is set. Number of converted points: " << number_of_converted_points << std::endl << std::endl;
			}

			void computing_visibility(
				CDT &cdt, 
				Log &log, 
				const Container_2D &input_2d, const size_t exec_step) {

				if (m_visibility.name() == "ray shooting" && m_pipeline_version == Pipeline_version::WITHOUT_SHAPE_DETECTION && !m_with_region_growing)
					assert(!"Ray shooting requires constrained edges!");

				if (m_visibility.name() == "blend") assert(!"Blend visibility is not worth trying!");

				// Compute visibility (0 - outside or 1 - inside) for each triangle in CDT above.
				std::cout << "(" << exec_step << ") visibility computation" << std::endl;

				m_visibility.save_info(m_visibility_save_info);
				m_visibility.set_approach(m_visibility_approach);
				m_visibility.set_method(m_visibility_method);
				m_visibility.set_number_of_samples(m_visibility_num_samples);
				m_visibility.show_progress(m_visibility_show_progress);
				m_visibility.set_norm_threshold(m_visibility_norm_threshold);
				m_visibility.set_number_of_neighbours(m_visibility_num_neighbours);
				m_visibility.set_sampler_type(m_visibility_sampler);
				m_visibility.set_number_of_rays_per_side(m_visibility_rays_per_side);
				m_visibility.set_small_edge_threshold(m_visibility_small_edge_threshold);
				m_visibility.set_angle_eps(m_visibility_angle_eps);

				const auto number_of_traversed_faces = m_visibility.compute(input_2d, cdt);
				log.out << "(" << exec_step << ") Visibility is computed. Number of traversed faces: " << number_of_traversed_faces << std::endl << std::endl;

				// Log eps_saver_wp; eps_saver_wp.save_visibility_eps(cdt, input, structured_points); // works only with basic test
				
				Log eps_saver; eps_saver.save_visibility_eps(cdt);
				Log ply_vis_saver; ply_vis_saver.save_cdt_ply(cdt, "tmp/visibility", "in");
			}

			void applying_graph_cut(
				CDT &cdt,
				Log &log, 
				const size_t exec_step) {

				// Apply graph cut.
				std::cout << "(" << exec_step << ") applying graph cut" << std::endl;

				m_graph_cut.save_info(m_graph_cut_save_info);
				m_graph_cut.set_alpha_parameter(m_graph_cut_alpha);
				m_graph_cut.set_beta_parameter(m_graph_cut_beta);
				m_graph_cut.set_gamma_parameter(m_graph_cut_gamma);

				m_graph_cut.max_flow(cdt);

				log.out << "(" << exec_step << ") Graph cut is applied." << std::endl << std::endl;
				Log ply_cdt_in; ply_cdt_in.save_cdt_ply(cdt, "tmp/after_cut", "in");
			}

			void splitting_buildings(
				Buildings &buildings, 
				CDT &cdt,
				Log &log, 
				const size_t exec_step) {

				// Split all buildings.
				std::cout << "(" << exec_step << ") splitting buildings" << std::endl;

				const auto number_of_buildings = m_building_splitter.split(cdt, buildings);

				log.out << "(" << exec_step << ") All buildings are found. Number of buildings: " << number_of_buildings << std::endl << std::endl;
			}

			void finding_buildings_boundaries(
				Buildings &buildings,
				Log &log, 
				const CDT &cdt, const size_t exec_step) {

				// Find building's walls.
				std::cout << "(" << exec_step << ") finding boundaries" << std::endl;

				m_building_outliner.save_info(m_building_boundaries_save_internal_info);
				m_building_outliner.set_max_inner_iterations(m_building_boundaries_max_inner_iters);
				m_building_outliner.set_max_outer_iterations(m_building_boundaries_max_outer_iters);
				m_building_outliner.set_boundary_type(m_building_boundary_type);
					
				m_building_outliner.find_boundaries(cdt, buildings);

				log.out << "(" << exec_step << ") All boundaries are found." << std::endl << std::endl; 

				// Log log_bounds; log_bounds.save_buildings_info(cdt, buildings, "tmp/buildings_info_with_boundaries");					
				// log_bounds.clear(); log_bounds.save_cdt_ply(cdt, "tmp/chosen_vertices"); // debugging info
			}

			void fitting_roofs(
				Buildings &buildings,
				Log &log, 
				const Plane_3 &fitted_ground_plane, const Face_points_map &fp_map, const Container_3D &input, const CDT &cdt, const size_t exec_step) {

				std::cout << "(" << exec_step << ") fitting roofs" << std::endl;
				switch (m_roof_fitter_type) {
					
					case Roof_fitter_type::MIN: 
						m_building_min_roof_fitter.fit_roof_heights(cdt, input, fp_map, fitted_ground_plane, buildings);
						break;

					case Roof_fitter_type::AVG: 
						m_building_avg_roof_fitter.fit_roof_heights(cdt, input, fp_map, fitted_ground_plane, buildings);
						break;

					case Roof_fitter_type::MAX: 
						m_building_max_roof_fitter.fit_roof_heights(cdt, input, fp_map, fitted_ground_plane, buildings);
						break;			

					default:
						assert(!"Wrong roof fitter type!");
						break;
				}
				log.out << "(" << exec_step << ") All roofs are fitted." << std::endl << std::endl;
				
				// Log log_roofs; log_roofs.save_buildings_info(cdt, buildings, "tmp/buildings_info_final");
			}

			void creating_lod0(
				Ground &ground_bbox, 
				Log &log, 
				const CDT &cdt, const Buildings &buildings, const Container_3D &input, const size_t exec_step) {

				// LOD0 reconstruction.

				m_lods.use_boundaries(m_use_boundaries);
				m_utils. template compute_ground_bbox<Ground, Ground_point>(input, ground_bbox);

				assert(!ground_bbox.empty());
				std::cout << "(" << exec_step << ") reconstructing lod0" << std::endl;

				Mesh mesh_0; Mesh_facet_colors mesh_facet_colors_0;
				m_lods.reconstruct_lod0(cdt, buildings, ground_bbox, mesh_0, mesh_facet_colors_0);

				log.out << "(" << exec_step << ") Final LOD0 is reconstructed." << std::endl << std::endl;
				Log lod_0_saver; lod_0_saver.save_mesh_as_ply(mesh_0, mesh_facet_colors_0, "LOD0");
			}

			void creating_lod1(
				Log &log, 
				const CDT &cdt, const Buildings &buildings, const Ground &ground_bbox, const size_t exec_step) {

				// LOD1 reconstruction.
				
				std::cout << "(" << exec_step << ") reconstructing lod1" << std::endl;

				Mesh mesh_1; Mesh_facet_colors mesh_facet_colors_1;
				m_lods.reconstruct_lod1(cdt, buildings, ground_bbox, mesh_1, mesh_facet_colors_1);

				log.out << "(" << exec_step << ") Final LOD1 is reconstructed." << std::endl;
				Log lod_1_saver; lod_1_saver.save_mesh_as_ply(mesh_1, mesh_facet_colors_1, "LOD1");
			}

			void finish_execution(
				Log &log, 
				const std::string &filename) {
				
				// Save log.
				std::cout << "... finishing\n" << std::endl;

				log.out << "\n\nFINISH EXECUTION";
				log.save(filename);
			}

		public:
			// Version 0.
			void create_lods_ver0() {

				// (--) ----------------------------------
				Log log; size_t exec_step = 0;
				start_execution(log);


				// (01) ----------------------------------
				Container_3D input;
				loading_data(input, log, ++exec_step);


				if (m_pipeline_version == Pipeline_version::WITH_SHAPE_DETECTION) {
					
					// (02) ---------------------------------- not used
					Planes all_planes;
					getting_all_planes(all_planes, log, input, ++exec_step);
				}


				// (03) ----------------------------------
				Indices ground_idxs, building_boundary_idxs, building_interior_idxs;
				applying_selection(ground_idxs, building_boundary_idxs, building_interior_idxs,
				log, input, ++exec_step);


				// (04) ----------------------------------
				Plane_3 base_ground_plane, fitted_ground_plane;
				ground_fitting(base_ground_plane, fitted_ground_plane, log, ground_idxs, input, ++exec_step);


				// (05) ----------------------------------					
				Boundary_data building_boundaries, boundary_clutter;					
				getting_boundary_points(building_boundaries, boundary_clutter, log, 
					building_boundary_idxs, building_interior_idxs, input, ++exec_step);


				if (m_pipeline_version == Pipeline_version::WITH_SHAPE_DETECTION) {

					// (06) ----------------------------------
					regularizing(building_boundaries, input, log, base_ground_plane, ++exec_step);
				}


				// (07) ----------------------------------
				Projected_points building_boundaries_projected, boundary_clutter_projected;
				projecting(building_boundaries_projected, boundary_clutter_projected, log, base_ground_plane, building_boundaries, boundary_clutter, input, ++exec_step);


				// (08) ----------------------------------
				if (m_pipeline_version == Pipeline_version::WITHOUT_SHAPE_DETECTION && m_with_region_growing) {

					if (m_use_grid_simplifier_first) {
						
						assert(m_clutter_new_point_type == Grid_new_point_type::CLOSEST);
						applying_grid_simplification(boundary_clutter_projected, log, ++exec_step);
					}

					detecting_2d_lines(boundary_clutter, boundary_clutter_projected, building_boundaries, building_boundaries_projected, log, input, ++exec_step);

					Lines lines; Segments segments;

					line_fitting(lines, log, building_boundaries, building_boundaries_projected, ++exec_step);
					creating_segments(segments, log, lines, building_boundaries, building_boundaries_projected, ++exec_step);
					applying_2d_structuring(log, lines, building_boundaries, building_boundaries_projected, ++exec_step);
				}


				if (m_pipeline_version == Pipeline_version::WITH_SHAPE_DETECTION) {

					// (09) ---------------------------------- not necessary in many cases
					if (m_clean_projected_points) 
						cleaning_projected_points(building_boundaries_projected, building_boundaries, log, ++exec_step);


					// (10) ----------------------------------
					Lines lines; 
					line_fitting(lines, log, building_boundaries, building_boundaries_projected, ++exec_step);


					// (11) ---------------------------------- not used
					Segments segments;
					creating_segments(segments, log, lines, building_boundaries, building_boundaries_projected, ++exec_step);


					// (12) ----------------------------------
					applying_2d_structuring(log, lines, building_boundaries, building_boundaries_projected, ++exec_step);
				}


				// (13) ----------------------------------
				if (m_pipeline_version == Pipeline_version::WITHOUT_SHAPE_DETECTION || 
				   (m_pipeline_version == Pipeline_version::WITH_SHAPE_DETECTION && m_add_cdt_clutter)) {

					if (!m_use_grid_simplifier_first)
						processing_clutter(boundary_clutter, boundary_clutter_projected, log, input, ++exec_step);
				}


				// (14) ----------------------------------
				CDT cdt;
				if (m_pipeline_version == Pipeline_version::WITHOUT_SHAPE_DETECTION && !m_with_region_growing) 
					m_add_cdt_clutter = true;
				creating_cdt(cdt, log, boundary_clutter, boundary_clutter_projected, input, ++exec_step);


				// (15) ----------------------------------
				Container_2D input_2d; Face_points_map fp_map;
				converting_3d_to_2d(input_2d, fp_map, log, cdt, input, ++exec_step);


				// (16) ----------------------------------
				computing_visibility(cdt, log, input_2d, ++exec_step);


				// (17) ----------------------------------
				applying_graph_cut(cdt, log, ++exec_step);

				
				// From now on we handle each building separately.

				// (18) ----------------------------------				
				Buildings buildings;
				splitting_buildings(buildings, cdt, log, ++exec_step);


				// (19) ----------------------------------				
				if (m_use_boundaries || m_building_boundary_type == Building_boundary_type::UNORIENTED) 
					finding_buildings_boundaries(buildings, log, cdt, ++exec_step);


				// (20) ----------------------------------
				fitting_roofs(buildings, log, fitted_ground_plane, fp_map, input, cdt, ++exec_step);


				// (21) ----------------------------------
				Ground ground_bbox;
				creating_lod0(ground_bbox, log, cdt, buildings, input, ++exec_step);


				// (22) ----------------------------------	
				creating_lod1(log, cdt, buildings, ground_bbox, ++exec_step);


				// (--) ----------------------------------
				finish_execution(log, "create_lods");
			}

		private:
			// Main components.
			Loader       m_loader;
			Preprocessor m_preprocessor;
			
			Building_boundary_selector m_building_boundary_selector;
			Building_interior_selector m_building_interior_selector;
			Clutter_selector           m_clutter_selector;
			Ground_selector            m_ground_selector;

			Vertical_regularizer m_vertical_regularizer;
			Ground_projector 	 m_ground_projector;
			Visibility_2 		 m_visibility;
			Utils 		 		 m_utils;
			Region_growing_2 	 m_region_growing;
			
			Graph_cut m_graph_cut;
			Lods m_lods;

			Building_splitter m_building_splitter;
			Building_outliner m_building_outliner;
			
			Building_min_roof_fitter m_building_min_roof_fitter;
			Building_avg_roof_fitter m_building_avg_roof_fitter;
			Building_max_roof_fitter m_building_max_roof_fitter;

			std::unique_ptr<Structuring_2> m_structuring;
			Clutter_processor m_clutter_processor;
			
			Grid_simplifier m_grid_simplifier;
			Thinning 		m_thinnning;


			// Global parameters.
			std::string m_prefix_path;
			std::string m_default_path;
			
			FT m_preprocessor_scale;
			FT m_structuring_epsilon;
			
			bool m_structuring_log;
			bool m_structuring_resample;
			bool m_structuring_get_all_points;

			bool m_add_cdt_clutter;
			bool m_add_cdt_bbox;

			bool m_visibility_save_info;
			Visibility_approach m_visibility_approach;
			Visibility_method m_visibility_method;
			size_t m_visibility_num_samples;

			bool m_graph_cut_save_info;
			FT m_graph_cut_alpha;
			FT m_graph_cut_beta;
			FT m_graph_cut_gamma;

			bool m_building_boundaries_save_internal_info;
			size_t m_building_boundaries_max_inner_iters;
			size_t m_building_boundaries_max_outer_iters;

			Roof_fitter_type m_roof_fitter_type;
			bool m_clean_projected_points;

			FT m_max_reg_angle;
			bool m_regularizer_reject_planes;
			bool m_use_boundaries;

			const Program_version  m_prog_version;
			
			Pipeline_version m_pipeline_version;

			bool m_visibility_show_progress;
			FT   m_visibility_norm_threshold;

			size_t 				   m_clutter_knn;
			FT 					   m_clutter_cell_length;
			Clutter_fitter_type    m_clutter_fitter_type;
			Clutter_new_point_type m_clutter_new_point_type;

			size_t m_visibility_num_neighbours;
			Visibility_sampler m_visibility_sampler;
			size_t m_visibility_rays_per_side;
			FT m_visibility_small_edge_threshold;

			Building_boundary_type m_building_boundary_type;
			FT m_visibility_angle_eps;

			Neighbour_search_type m_thinning_neighbour_search_type;
			FT 				      m_thinning_fuzzy_radius;
			Thinning_type 		  m_thinning_type;
			
			FT 	   m_region_growing_epsilon;
			FT     m_region_growing_cluster_epsilon;
			FT 	   m_region_growing_normal_threshold;
			size_t m_region_growing_min_points;

			bool m_with_region_growing;
			bool m_use_grid_simplifier_first;

			FT m_alpha_shape_size;
			bool m_use_alpha_shapes;

			Structuring_corner_algorithm 		   m_structuring_corner_algorithm;
			Structuring_adjacency_threshold_method m_structuring_adjacency_method;
			FT 									   m_structuring_adjacency_value;
			bool m_structuring_global_everywhere;


			// Assert default values of all global parameters.
			void assert_global_parameters() {
			
				assert(m_default_path != "default");

				assert(m_preprocessor_scale  != -FT(1));
				assert(m_structuring_epsilon != -FT(1));

				assert(!m_add_cdt_bbox);
				assert(m_visibility_num_samples != 0);

				assert(!(m_visibility_approach == Visibility_approach::FACE_BASED  && m_visibility_method == Visibility_method::POINT_BASED_CLASSIFICATION));
				assert(!(m_visibility_approach == Visibility_approach::POINT_BASED && m_visibility_method == Visibility_method::FACE_BASED_COUNT));
				assert(!(m_visibility_approach == Visibility_approach::POINT_BASED && m_visibility_method == Visibility_method::FACE_BASED_NATURAL_NEIGHBOURS));
				assert(!(m_visibility_approach == Visibility_approach::POINT_BASED && m_visibility_method == Visibility_method::FACE_BASED_BARYCENTRIC));

				assert(m_graph_cut_alpha != -FT(1));
				assert(m_graph_cut_beta  != -FT(1));
				assert(m_graph_cut_gamma != -FT(1));

				assert(m_building_boundaries_max_inner_iters != 0);
				assert(m_building_boundaries_max_outer_iters != 0);

				assert(m_max_reg_angle != -FT(1));
				assert(m_visibility_norm_threshold != -FT(1));
				
				assert(m_clutter_knn > 1);
				assert(m_clutter_cell_length != -FT(1));

				assert(m_visibility_num_neighbours > 1);
				assert(m_visibility_rays_per_side  > 0);
				assert(m_visibility_angle_eps != -FT(1));

				assert(m_thinning_fuzzy_radius != -FT(1));
				assert(m_structuring_adjacency_value > FT(0));

				if (m_with_region_growing) {
					assert(m_region_growing_epsilon 		 != -FT(1));
					assert(m_region_growing_cluster_epsilon  != -FT(1));
					assert(m_region_growing_normal_threshold != -FT(1));
					assert(m_region_growing_min_points 		 !=     0);
				}

				if (m_use_alpha_shapes) {
					assert(m_alpha_shape_size > FT(0));
				}

				assert(m_prefix_path != "path_to_the_data_folder");
			}


			// Set all global parameters.
			void set_global_parameters() {

				// General parameters. Not important!
				m_add_cdt_bbox = false;
				
				m_structuring_log 	   	 				 = false;
				m_visibility_save_info 					 = false;
				m_graph_cut_save_info 					 = false;
				m_building_boundaries_save_internal_info = false;

				m_building_boundaries_max_inner_iters = 1000;
				m_building_boundaries_max_outer_iters = 1000000;

				m_visibility_show_progress  = true;
				m_visibility_norm_threshold = 1000.0;


				// More important.
				m_regularizer_reject_planes  = true;
				m_structuring_resample 	 	 = true;
				m_structuring_get_all_points = false;
				m_clean_projected_points 	 = true;
				
				m_roof_fitter_type 			   = Roof_fitter_type::AVG;
				m_structuring_corner_algorithm = Structuring_corner_algorithm::GRAPH_BASED;
				m_structuring_adjacency_method = Structuring_adjacency_threshold_method::LOCAL;
				m_structuring_adjacency_value  = 0.00001;

				m_preprocessor_scale  		= 2.0;
				m_clutter_fitter_type 		= Clutter_fitter_type::LINE;	
				m_use_grid_simplifier_first = false;
				
				m_visibility_num_neighbours 	  = 6;
				m_visibility_rays_per_side  	  = 10;
				m_visibility_small_edge_threshold = -1000000.0; // not used
				m_structuring_global_everywhere   = true;

				m_graph_cut_alpha = 1.0;    // should not change anything but should be bigger or equal to 1
				m_graph_cut_gamma = 1000.0; // is not used in the pipeline without shape detection (or without structuring), otherwise should be some big value


				// The most important!
				const Main_test_data_type test_data_type = Main_test_data_type::PARIS_ETH;
				switch (test_data_type) {

					case Main_test_data_type::BASIC:
						set_basic_parameters();
						break;

					case Main_test_data_type::COMPLEX:
						set_complex_parameters();
						break;

					case Main_test_data_type::P10:
						set_p10_parameters();
						break;

					case Main_test_data_type::PARIS:
						set_paris_parameters();
						break;

					case Main_test_data_type::PARIS_FULL:
						set_paris_full_parameters();
						break;

					case Main_test_data_type::PARIS_ETH:
						set_paris_eth_parameters();
						break;

					case Main_test_data_type::PARIS_FULL_ETH:
						set_paris_full_eth_parameters();
						break;

					case Main_test_data_type::RESIDENT_TILE_1:
						set_resident_tile_1_parameters();
						break;

					case Main_test_data_type::RESIDENT_TILE_2:
						set_resident_tile_2_parameters();
						break;

					case Main_test_data_type::RESIDENT_TILE_3:
						set_resident_tile_3_parameters();
						break;

					case Main_test_data_type::PARIS_TILE_1:
						set_paris_tile_1_parameters();
						break;

					case Main_test_data_type::PARIS_TILE_2:
						set_paris_tile_2_parameters();
						break;

					default:
						assert(!"Wrong test data!");
						break;
				}
			}


			// no shape detection
			void set_basic_parameters() {

				// To load basic parameters from stub, add stub to the loader class in LOD_traits!
				// Stub works with these basic parameters. Actually it is the same data set.

				// All main parameters are set below.
				m_default_path     = m_prefix_path + "basic_test/data";
				m_pipeline_version = Pipeline_version::WITH_SHAPE_DETECTION;

				m_visibility_approach 	 		 = Visibility_approach::POINT_BASED;
				m_visibility_method   	 		 = Visibility_method::POINT_BASED_CLASSIFICATION; // use natural neighbours for without_shape_detection
				m_visibility_sampler 	 		 = Visibility_sampler::BARYCENTRE;
				m_thinning_neighbour_search_type = Neighbour_search_type::KNN;
				m_building_boundary_type 		 = Building_boundary_type::ORIENTED;
				m_clutter_new_point_type 		 = Clutter_new_point_type::BARYCENTRE; // BARYCENTRE - keeps average position of the removed points, CENTROID - inserts new point in the centre of the grid cell, CLOSEST - to the barycentre	
				m_thinning_type 	  			 = Thinning_type::NAIVE;

				m_thinning_fuzzy_radius  = 0.000001; // if zero, is not used; the bigger radius, the longer it works
				m_visibility_angle_eps   = 0.0; 	 // used for removing thin triangles; if zero, is not used
				m_max_reg_angle          = 10.0;	 // in average should be 10-20 degrees.
				m_structuring_epsilon 	 = 0.025; 	 // the most important parameter!!! Depends on the dataset.
				m_add_cdt_clutter     	 = true;	 // is always true if shape detection is not used
				m_visibility_num_samples = 1;		 // the more samples, the slower but better quality in visibility
				m_graph_cut_beta 		 = 100000.0; // smaller value for less inside triangles
				m_clutter_knn 			 = 2;		 // the smaller value, the less thinning is performed
				m_clutter_cell_length    = 0.025;	 // the bigger value, the more points are removed in the grid simplify
				m_use_boundaries 		 = true;     // use or not outliner to build walls
				
				m_with_region_growing 	 		  = false; // use 2D region growing for structuring clutter
				m_region_growing_epsilon 		  = 0.0;   // distance to the line
				m_region_growing_cluster_epsilon  = 0.0;   // distance between neighbouring points
				m_region_growing_normal_threshold = 0.0;   // difference between the line normal and the point normal
				m_region_growing_min_points 	  = 0;     // min number of points per shape

				m_use_alpha_shapes = false; // if true, we use alpha shapes to extract boundary points from building roofs
				m_alpha_shape_size = -1.0;  // size of the ball used in alpha shapes to extract boundaries
			}


			// open cv random forest
			void set_complex_parameters() {

				// SWITCH TO RAY SHOOTING HERE!

				// All main parameters are set below.
				// If using ray shooting here, we need to use with_shape_detection.
				m_default_path     = m_prefix_path + "complex_test/data_region_growing";
				m_pipeline_version = Pipeline_version::WITH_SHAPE_DETECTION;

				m_visibility_approach 	 		 = Visibility_approach::POINT_BASED;
				m_visibility_method   	 		 = Visibility_method::POINT_BASED_CLASSIFICATION;
				m_visibility_sampler 	 		 = Visibility_sampler::UNIFORM_SUBDIVISION;
				m_thinning_neighbour_search_type = Neighbour_search_type::KNN;
				m_building_boundary_type 		 = Building_boundary_type::ORIENTED;
				m_clutter_new_point_type 		 = Clutter_new_point_type::BARYCENTRE;
				m_thinning_type 	  			 = Thinning_type::NAIVE;

				m_thinning_fuzzy_radius  = 0.001;
				m_visibility_angle_eps   = 0.0; 	 
				m_max_reg_angle          = 10.0;
				m_structuring_epsilon 	 = 0.0005;
				m_add_cdt_clutter     	 = false;
				m_visibility_num_samples = 3;
				m_graph_cut_beta 		 = 100000.0; // use 1.0 for without_shape_detection
				m_clutter_knn 			 = 12;
				m_clutter_cell_length    = 0.015;
				m_use_boundaries 		 = true;
				
				m_with_region_growing 	 		  = false;
				m_region_growing_epsilon 		  = 0.0;  
				m_region_growing_cluster_epsilon  = 0.0;  
				m_region_growing_normal_threshold = 0.0;  
				m_region_growing_min_points 	  = 0;  

				m_use_alpha_shapes = false;
				m_alpha_shape_size = -1.0;   
			}


			// weighted sum
			void set_p10_parameters() {

				// YOU CAN USE HERE RAY SHOOTING FOR WITH_SHAPE_DETECTION!

				// All main parameters are set below.
				m_default_path     = m_prefix_path + "p10_test/data_region_growing_weighted_sum";
				m_pipeline_version = Pipeline_version::WITH_SHAPE_DETECTION;

				m_visibility_approach 	 		 = Visibility_approach::FACE_BASED;
				m_visibility_method   	 		 = Visibility_method::FACE_BASED_NATURAL_NEIGHBOURS;
				m_visibility_sampler 	 		 = Visibility_sampler::UNIFORM_SUBDIVISION;
				m_thinning_neighbour_search_type = Neighbour_search_type::CIRCLE;
				m_building_boundary_type 		 = Building_boundary_type::ORIENTED;
				m_clutter_new_point_type 		 = Clutter_new_point_type::BARYCENTRE;
				m_thinning_type 	  			 = Thinning_type::NAIVE;

				m_thinning_fuzzy_radius  = 5.0;
				m_visibility_angle_eps   = 0.0; 
				m_max_reg_angle          = 15.0;
				m_structuring_epsilon 	 = 0.2;
				m_add_cdt_clutter     	 = false;
				m_visibility_num_samples = 2;
				m_graph_cut_beta 		 = 100000.0; // 10.0 for without_shape_detection
				m_clutter_knn 			 = 12;
				m_clutter_cell_length    = 4.0;
				m_use_boundaries 		 = true;
				
				m_with_region_growing 	 		  = false;
				m_region_growing_epsilon 		  = 0.0;  
				m_region_growing_cluster_epsilon  = 0.0;  
				m_region_growing_normal_threshold = 0.0;  
				m_region_growing_min_points 	  = 0;   

				m_use_alpha_shapes = false;
				m_alpha_shape_size = -1.0; 
			}


			// weighted sum
			void set_paris_parameters() {

				// All main parameters are set below.
				m_default_path     = m_prefix_path + "paris_test/data_region_growing_weighted_sum";
				m_pipeline_version = Pipeline_version::WITHOUT_SHAPE_DETECTION;

				m_visibility_approach 	 		 = Visibility_approach::FACE_BASED;
				m_visibility_method   	 		 = Visibility_method::FACE_BASED_NATURAL_NEIGHBOURS; // point based for with_shape_detection
				m_visibility_sampler 	 		 = Visibility_sampler::UNIFORM_SUBDIVISION;
				m_thinning_neighbour_search_type = Neighbour_search_type::CIRCLE;
				m_building_boundary_type 		 = Building_boundary_type::ORIENTED;
				m_clutter_new_point_type 		 = Clutter_new_point_type::BARYCENTRE;
				m_thinning_type 	  			 = Thinning_type::NAIVE;

				m_thinning_fuzzy_radius  = 5.0;
				m_visibility_angle_eps   = 0.001; 
				m_max_reg_angle          = 10.0;
				m_structuring_epsilon 	 = 1.5;
				m_add_cdt_clutter     	 = true;
				m_visibility_num_samples = 1;
				m_graph_cut_beta 		 = 35.0; // 15.0 for with_shape_detection
				m_clutter_knn 			 = 12;
				m_clutter_cell_length    = 10.0;
				m_use_boundaries 		 = true;
				
				m_with_region_growing 	 		  = false;
				m_region_growing_epsilon 		  = 0.0;  
				m_region_growing_cluster_epsilon  = 0.0;  
				m_region_growing_normal_threshold = 0.0;  
				m_region_growing_min_points 	  = 0; 

				m_use_alpha_shapes = false;
				m_alpha_shape_size = -1.0;    
			}


			// weighted sum
			void set_paris_full_parameters() {

				// All main parameters are set below.
				m_default_path     = m_prefix_path + "paris_full_test/data_region_growing_weighted_sum";
				m_pipeline_version = Pipeline_version::WITHOUT_SHAPE_DETECTION;

				m_visibility_approach 	 		 = Visibility_approach::FACE_BASED;
				m_visibility_method   	 		 = Visibility_method::FACE_BASED_NATURAL_NEIGHBOURS; // point based for with_shape_detection
				m_visibility_sampler 	 		 = Visibility_sampler::UNIFORM_SUBDIVISION;
				m_thinning_neighbour_search_type = Neighbour_search_type::CIRCLE;
				m_building_boundary_type 		 = Building_boundary_type::UNORIENTED;
				m_clutter_new_point_type 		 = Clutter_new_point_type::BARYCENTRE;
				m_thinning_type 	  			 = Thinning_type::NAIVE;

				m_thinning_fuzzy_radius  = 5.0;
				m_visibility_angle_eps   = 0.001; 
				m_max_reg_angle          = 10.0;
				m_structuring_epsilon 	 = 1.5;
				m_add_cdt_clutter     	 = true;
				m_visibility_num_samples = 1;
				m_graph_cut_beta 		 = 35.0; // 15.0 for with_shape_detection
				m_clutter_knn 			 = 12;
				m_clutter_cell_length    = 10.0;
				m_use_boundaries 		 = true;
				
				m_with_region_growing 	 		  = false;
				m_region_growing_epsilon 		  = 0.0;  
				m_region_growing_cluster_epsilon  = 0.0;  
				m_region_growing_normal_threshold = 0.0;  
				m_region_growing_min_points 	  = 0;   

				m_use_alpha_shapes = false;
				m_alpha_shape_size = -1.0; 
			}


			// eth random forest
			void set_paris_eth_parameters() {

				// All main parameters are set below.
				m_default_path     = m_prefix_path + "paris_test/data_region_growing_eth";
				m_pipeline_version = Pipeline_version::WITHOUT_SHAPE_DETECTION;

				m_visibility_approach 	 		 = Visibility_approach::FACE_BASED;
				m_visibility_method   	 		 = Visibility_method::FACE_BASED_NATURAL_NEIGHBOURS; // point based for with_shape_detection
				m_visibility_sampler 	 		 = Visibility_sampler::UNIFORM_SUBDIVISION;
				m_thinning_neighbour_search_type = Neighbour_search_type::CIRCLE;
				m_building_boundary_type 		 = Building_boundary_type::UNORIENTED;
				m_clutter_new_point_type 		 = Clutter_new_point_type::CLOSEST;
				m_thinning_type 	  			 = Thinning_type::NAIVE;
				m_structuring_adjacency_method 	 = Structuring_adjacency_threshold_method::GLOBAL;

				m_thinning_fuzzy_radius  = 5.0;
				m_visibility_angle_eps   = 0.18; 
				m_max_reg_angle          = 10.0;
				m_structuring_epsilon 	 = 5.0;
				m_add_cdt_clutter     	 = false;
				m_visibility_num_samples = 1;
				m_graph_cut_beta 		 = 100000.0; // 35.0 with_clutter // 15.0 for with_shape_detection
				m_clutter_knn 			 = 12;
				m_clutter_cell_length    = 1.3;
				m_use_boundaries 		 = true;

				m_use_grid_simplifier_first = true;
				m_with_region_growing 	 	= true;

				m_region_growing_epsilon 		  = 3.2;
				m_region_growing_cluster_epsilon  = 2.9;
				m_region_growing_normal_threshold = 0.7;  
				m_region_growing_min_points 	  = 10;

				m_use_alpha_shapes = true;
				m_alpha_shape_size = 5.0;
				m_graph_cut_gamma  = 10000.0;

				m_structuring_get_all_points  = true;
				m_structuring_adjacency_value = 12.0;
			}


			// eth random forest
			void set_paris_full_eth_parameters() {

				// All main parameters are set below.
				m_default_path     = m_prefix_path + "paris_full_test/data_region_growing_eth";
				m_pipeline_version = Pipeline_version::WITHOUT_SHAPE_DETECTION;

				m_visibility_approach 	 		 = Visibility_approach::FACE_BASED;
				m_visibility_method   	 		 = Visibility_method::FACE_BASED_NATURAL_NEIGHBOURS;
				m_visibility_sampler 	 		 = Visibility_sampler::UNIFORM_SUBDIVISION;
				m_thinning_neighbour_search_type = Neighbour_search_type::CIRCLE;
				m_building_boundary_type 		 = Building_boundary_type::UNORIENTED;
				m_clutter_new_point_type 		 = Clutter_new_point_type::CLOSEST;
				m_thinning_type 	  			 = Thinning_type::NAIVE;
				m_structuring_adjacency_method 	 = Structuring_adjacency_threshold_method::GLOBAL;
				m_structuring_corner_algorithm   = Structuring_corner_algorithm::NO_T_CORNERS;

				m_thinning_fuzzy_radius  = 5.0;
				m_visibility_angle_eps   = 0.18; 
				m_max_reg_angle          = 10.0;
				m_structuring_epsilon 	 = 5.0;
				m_add_cdt_clutter     	 = false;
				m_visibility_num_samples = 1;
				m_graph_cut_beta 		 = 100000.0;
				m_clutter_knn 			 = 12;
				m_clutter_cell_length    = 1.3;
				m_use_boundaries 		 = true;

				m_use_grid_simplifier_first = true;
				m_with_region_growing 	 	= true;

				m_region_growing_epsilon 		  = 3.2;
				m_region_growing_cluster_epsilon  = 2.9;
				m_region_growing_normal_threshold = 0.7;  
				m_region_growing_min_points 	  = 10;

				m_use_alpha_shapes = true;
				m_alpha_shape_size = 5.0;

				m_structuring_get_all_points  = true;
				m_structuring_adjacency_value = 12.0;
			}


			// most likely it is not going to look good for the moment
			void set_resident_tile_1_parameters() {

				// All main parameters are set below.
				m_default_path     = m_prefix_path + "residential_test/tile_1/data_region_growing_eth";
				m_pipeline_version = Pipeline_version::WITHOUT_SHAPE_DETECTION;

				m_visibility_approach 	 		 = Visibility_approach::POINT_BASED;
				m_visibility_method   	 		 = Visibility_method::POINT_BASED_CLASSIFICATION;
				m_visibility_sampler 	 		 = Visibility_sampler::UNIFORM_SUBDIVISION;
				m_thinning_neighbour_search_type = Neighbour_search_type::CIRCLE;
				m_building_boundary_type 		 = Building_boundary_type::UNORIENTED;
				m_clutter_new_point_type 		 = Clutter_new_point_type::CLOSEST;
				m_thinning_type 	  			 = Thinning_type::NAIVE;
				m_structuring_adjacency_method 	 = Structuring_adjacency_threshold_method::GLOBAL;
				m_structuring_corner_algorithm   = Structuring_corner_algorithm::NO_T_CORNERS;

				m_thinning_fuzzy_radius  = 5.0;
				m_visibility_angle_eps   = 0.18;
				m_max_reg_angle          = 10.0;
				m_structuring_epsilon 	 = 2.0;
				m_add_cdt_clutter     	 = false;
				m_visibility_num_samples = 1;
				m_graph_cut_beta 		 = 100000.0;
				m_clutter_knn 			 = 12;
				m_clutter_cell_length    = 1.3;
				m_use_boundaries 		 = true;

				m_use_grid_simplifier_first = true;
				m_with_region_growing 	 	= true;

				m_region_growing_epsilon 		  = 3.2;
				m_region_growing_cluster_epsilon  = 2.9;
				m_region_growing_normal_threshold = 0.7;  
				m_region_growing_min_points 	  = 10;

				m_use_alpha_shapes = true;
				m_alpha_shape_size = 5.0;

				m_structuring_get_all_points    = true;
				m_structuring_adjacency_value   = 4.0;
				m_structuring_global_everywhere = false;
			}

			// does not look good for the moment
			void set_resident_tile_2_parameters() {

				// All main parameters are set below.
				m_default_path     = m_prefix_path + "residential_test/tile_2/data_region_growing_eth";
				m_pipeline_version = Pipeline_version::WITHOUT_SHAPE_DETECTION;

				m_visibility_approach 	 		 = Visibility_approach::FACE_BASED;
				m_visibility_method   	 		 = Visibility_method::FACE_BASED_NATURAL_NEIGHBOURS;
				m_visibility_sampler 	 		 = Visibility_sampler::UNIFORM_SUBDIVISION;
				m_thinning_neighbour_search_type = Neighbour_search_type::CIRCLE;
				m_building_boundary_type 		 = Building_boundary_type::UNORIENTED;
				m_clutter_new_point_type 		 = Clutter_new_point_type::CLOSEST;
				m_thinning_type 	  			 = Thinning_type::NAIVE;
				m_structuring_adjacency_method 	 = Structuring_adjacency_threshold_method::GLOBAL;
				m_structuring_corner_algorithm   = Structuring_corner_algorithm::NO_T_CORNERS;

				m_thinning_fuzzy_radius  = 5.0;
				m_visibility_angle_eps   = 0.18;
				m_max_reg_angle          = 10.0;
				m_structuring_epsilon 	 = 2.0;
				m_add_cdt_clutter     	 = false;
				m_visibility_num_samples = 1;
				m_graph_cut_beta 		 = 100000.0;
				m_clutter_knn 			 = 12;
				m_clutter_cell_length    = 1.3;
				m_use_boundaries 		 = true;

				m_use_grid_simplifier_first = true;
				m_with_region_growing 	 	= true;

				m_region_growing_epsilon 		  = 3.2;
				m_region_growing_cluster_epsilon  = 2.9;
				m_region_growing_normal_threshold = 0.7;  
				m_region_growing_min_points 	  = 10;

				m_use_alpha_shapes = true;
				m_alpha_shape_size = 5.0;

				m_structuring_get_all_points    = true;
				m_structuring_adjacency_value   = 4.0;
				m_structuring_global_everywhere = false;
			}


			// residential area - many small buildings and two large buildings, one tile, eth random forest
			void set_resident_tile_3_parameters() {

				// All main parameters are set below.
				m_default_path     = m_prefix_path + "residential_test/tile_3/data_region_growing_eth";
				m_pipeline_version = Pipeline_version::WITHOUT_SHAPE_DETECTION;

				m_visibility_approach 	 		 = Visibility_approach::FACE_BASED;
				m_visibility_method   	 		 = Visibility_method::FACE_BASED_NATURAL_NEIGHBOURS;
				m_visibility_sampler 	 		 = Visibility_sampler::UNIFORM_SUBDIVISION;
				m_thinning_neighbour_search_type = Neighbour_search_type::CIRCLE;
				m_building_boundary_type 		 = Building_boundary_type::UNORIENTED;
				m_clutter_new_point_type 		 = Clutter_new_point_type::CLOSEST;
				m_thinning_type 	  			 = Thinning_type::NAIVE;
				m_structuring_adjacency_method 	 = Structuring_adjacency_threshold_method::GLOBAL;

				m_thinning_fuzzy_radius  = 5.0;
				m_visibility_angle_eps   = 0.18;
				m_max_reg_angle          = 10.0;
				m_structuring_epsilon 	 = 2.0;
				m_add_cdt_clutter     	 = false;
				m_visibility_num_samples = 1;
				m_graph_cut_beta 		 = 100000.0;
				m_clutter_knn 			 = 12;
				m_clutter_cell_length    = 1.3;
				m_use_boundaries 		 = true;

				m_use_grid_simplifier_first = true;
				m_with_region_growing 	 	= true;

				m_region_growing_epsilon 		  = 3.2;
				m_region_growing_cluster_epsilon  = 2.9;
				m_region_growing_normal_threshold = 0.7;  
				m_region_growing_min_points 	  = 3;

				m_use_alpha_shapes = true;
				m_alpha_shape_size = 5.0;

				m_structuring_get_all_points    = true;
				m_structuring_adjacency_value   = 4.0;
				m_structuring_global_everywhere = false;
			}


			// eth random forest
			void set_paris_tile_1_parameters() {

				// All main parameters are set below.
				m_default_path     = m_prefix_path + "paris_tiles_test/tile_1/data_region_growing_eth";
				m_pipeline_version = Pipeline_version::WITHOUT_SHAPE_DETECTION;

				m_visibility_approach 	 		 = Visibility_approach::FACE_BASED;
				m_visibility_method   	 		 = Visibility_method::FACE_BASED_NATURAL_NEIGHBOURS;
				m_visibility_sampler 	 		 = Visibility_sampler::UNIFORM_SUBDIVISION;
				m_thinning_neighbour_search_type = Neighbour_search_type::CIRCLE;
				m_building_boundary_type 		 = Building_boundary_type::UNORIENTED;
				m_clutter_new_point_type 		 = Clutter_new_point_type::CLOSEST;
				m_thinning_type 	  			 = Thinning_type::NAIVE;
				m_structuring_adjacency_method 	 = Structuring_adjacency_threshold_method::GLOBAL;
				m_structuring_corner_algorithm   = Structuring_corner_algorithm::NO_T_CORNERS;

				m_thinning_fuzzy_radius  = 5.0;
				m_visibility_angle_eps   = 0.18; 
				m_max_reg_angle          = 10.0;
				m_structuring_epsilon 	 = 5.0;
				m_add_cdt_clutter     	 = false;
				m_visibility_num_samples = 1;
				m_graph_cut_beta 		 = 100000.0;
				m_clutter_knn 			 = 12;
				m_clutter_cell_length    = 1.3;
				m_use_boundaries 		 = true;

				m_use_grid_simplifier_first = true;
				m_with_region_growing 	 	= true;

				m_region_growing_epsilon 		  = 3.2;
				m_region_growing_cluster_epsilon  = 2.9;
				m_region_growing_normal_threshold = 0.7;  
				m_region_growing_min_points 	  = 10;

				m_use_alpha_shapes = true;
				m_alpha_shape_size = 5.0;

				m_structuring_get_all_points  = true;
				m_structuring_adjacency_value = 12.0;
			}


			// eth random forest
			void set_paris_tile_2_parameters() {

				// All main parameters are set below.
				m_default_path     = m_prefix_path + "paris_tiles_test/tile_2/data_region_growing_eth";
				m_pipeline_version = Pipeline_version::WITHOUT_SHAPE_DETECTION;

				m_visibility_approach 	 		 = Visibility_approach::FACE_BASED;
				m_visibility_method   	 		 = Visibility_method::FACE_BASED_NATURAL_NEIGHBOURS;
				m_visibility_sampler 	 		 = Visibility_sampler::UNIFORM_SUBDIVISION;
				m_thinning_neighbour_search_type = Neighbour_search_type::CIRCLE;
				m_building_boundary_type 		 = Building_boundary_type::UNORIENTED;
				m_clutter_new_point_type 		 = Clutter_new_point_type::CLOSEST;
				m_thinning_type 	  			 = Thinning_type::NAIVE;
				m_structuring_adjacency_method 	 = Structuring_adjacency_threshold_method::GLOBAL;
				m_structuring_corner_algorithm   = Structuring_corner_algorithm::NO_T_CORNERS;

				m_thinning_fuzzy_radius  = 5.0;
				m_visibility_angle_eps   = 0.18; 
				m_max_reg_angle          = 10.0;
				m_structuring_epsilon 	 = 5.0;
				m_add_cdt_clutter     	 = false;
				m_visibility_num_samples = 1;
				m_graph_cut_beta 		 = 100000.0;
				m_clutter_knn 			 = 12;
				m_clutter_cell_length    = 1.3;
				m_use_boundaries 		 = true;

				m_use_grid_simplifier_first = true;
				m_with_region_growing 	 	= true;

				m_region_growing_epsilon 		  = 3.2;
				m_region_growing_cluster_epsilon  = 2.9;
				m_region_growing_normal_threshold = 0.7;  
				m_region_growing_min_points 	  = 10;

				m_use_alpha_shapes = true;
				m_alpha_shape_size = 5.0;

				m_structuring_get_all_points  = true;
				m_structuring_adjacency_value = 12.0;
			}


			//////////////////////
			// Not used functions!

			void projecting_without_shape_detection(
				Projected_points &building_boundaries_projected, 
				Log &log, 
				const Plane_3 &base_ground_plane, const Indices &building_boundary_idxs, const Container_3D &input, const size_t exec_step) {
				
				// Project all boundary points onto the ground plane.
				std::cout << "(" << exec_step << ") projecting" << std::endl;

				const auto number_of_projected_points = m_ground_projector.project_with_indices(input, building_boundary_idxs, base_ground_plane, building_boundaries_projected);
				
				Log proj_saver; proj_saver.export_projected_points_as_xyz("tmp/projected_boundaries", building_boundaries_projected, m_default_path);
				log.out << "(" << exec_step << ") Building's boundary points are projected. Number of projected points: " << number_of_projected_points << std::endl << std::endl;
			}

			void creating_cdt_without_shape_detection(
				CDT &cdt, 
				Log &log, 
				const Projected_points &building_boundaries_projected, const size_t exec_step) {

				// Compute constrained Delaunay triangulation of the input boundary points.
				std::cout << "(" << exec_step << ") creating cdt" << std::endl;
 
				const auto number_of_faces = m_utils.compute_delaunay(building_boundaries_projected, cdt);

				log.out << "(" << exec_step << ") Delaunay triangulation of the input boundary points is built. Number of faces: " << number_of_faces << std::endl << std::endl;
			}
		};
	}
}

#endif // CGAL_LEVEL_OF_DETAIL_BASE_H	