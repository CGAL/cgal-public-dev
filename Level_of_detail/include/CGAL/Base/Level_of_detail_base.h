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

			typedef typename Traits::Structuring_2 	Structuring_2;
			typedef typename Traits::Visibility_2  	Visibility_2;
			
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


			//////////////
			// Main class!
			Level_of_detail_base() :
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
			m_reject_planes(true), 
			m_use_boundaries(true)
			{ } // Do I need to create an instance of these traits here?


			//////////////////
			// Main functions!

			// Version a.
			void create_lods_a() {

				// (START) Create log.
				std::cout << "\nstarting ..." << std::endl;
				Log log; log.out << "START EXECUTION\n\n\n";

				set_global_parameters();
				assert_global_parameters();


				// ----------------------------------

				// (1) Read data.
				std::cout << "(1) loading" << std::endl;

				Container_3D input;
				m_loader.get_data(m_default_path + ".ply", input);

				log.out << "(1) Data are loaded. Number of points: " << input.number_of_points() << std::endl << std::endl;

				// Log mock_saver; mock_saver.save_ply<Traits, Container_3D>(input, "basic_mock", true);


				// ----------------------------------

				// (2) Find a set of planes related to the points. Basically here we emulate RANSAC.
				// For each plane we store indices of all points contained in this plane.
				std::cout << "(2) planes" << std::endl;

				Planes all_planes;
				const auto number_of_planes = m_preprocessor.get_planes(input, all_planes);
				assert(number_of_planes >= 0);

				log.out << "(2) Planes are found. Number of planes: " << number_of_planes << std::endl << std::endl;


				// ----------------------------------

				// (3) Split data with respect to 4 different semantic labels.
				std::cout << "(3) selection" << std::endl;

				Indices building_boundary_idxs, building_interior_idxs, clutter_idxs, ground_idxs;

				m_clutter_selector.select_elements(input, std::back_inserter(clutter_idxs));
				m_ground_selector.select_elements( input, std::back_inserter(ground_idxs));

				m_building_boundary_selector.select_elements(input, std::back_inserter(building_boundary_idxs));
				m_building_interior_selector.select_elements(input, std::back_inserter(building_interior_idxs));

				log.out << "(3) Clutter is found. Number of elements: " 			 << clutter_idxs.size() << std::endl;
				log.out << "(.) Ground is found. Number of elements: " 			     << ground_idxs.size() << std::endl;
				log.out << "(.) Building boundaries are found. Number of elements: " << building_boundary_idxs.size() << std::endl;
				log.out << "(3) Building interiors are found. Number of elements: "  << building_interior_idxs.size() << std::endl << std::endl;


				// ----------------------------------

				// (4) Create plane from the ground points.
				std::cout << "(4) ground plane fitting" << std::endl;

				Plane_3 ground_fitted_plane;
				m_utils.fit_ground_plane(input, ground_idxs, ground_fitted_plane);

				Plane_3 ground_plane = Plane_3(FT(0), FT(0), FT(1), FT(0)); // use XY plane instead
				log.out << "(4) Ground plane is fitted: " << ground_plane << std::endl << std::endl;


				// ----------------------------------

				// (5) Map indices from all detected planes to the ones that are a part of the given facades.
				std::cout << "(5) getting boundaries" << std::endl;

				Boundary_data building_boundaries, boundary_clutter;
				const auto number_of_boundaries = m_preprocessor.get_boundaries(input, building_boundary_idxs, building_boundaries, boundary_clutter);

				log.out << "(5) Planes for building's boundary are found. Number of planes: " << number_of_boundaries << std::endl;
				log.out << "(5) Boundary clutter is found. Number of points: " << boundary_clutter.at(0).size() << std::endl << std::endl;


				// ----------------------------------

				// (6) Make all nearly vertical planes in the building's boundary exactly vertical.
				std::cout << "(6) regularizing" << std::endl;

				m_vertical_regularizer.set_max_regularization_angle(m_max_reg_angle);
				m_vertical_regularizer.reject_planes(m_reject_planes);

				const auto number_of_regularized_planes = m_vertical_regularizer.regularize(building_boundaries, input, ground_plane);

				log.out << "(6) Building's nearly vertical planes are regularized. Number of regularized planes: " << number_of_regularized_planes <<
				", number of rejected planes: " << number_of_boundaries - building_boundaries.size() << std::endl << std::endl;

				// Log ply_saver; ply_saver.save_ply<Kernel, Container_3D>(input, "regularized", true);


				// ----------------------------------

				// (7) Project all vertical building's boundaries onto the ground plane.
				std::cout << "(7) projecting" << std::endl;

				Projected_points building_boundaries_projected; 
				auto number_of_projected_points = m_ground_projector.project_with_planes(input, building_boundaries, ground_plane, building_boundaries_projected);
				log.out << "(7) Building's boundary planar points are projected. Number of projected points: " << number_of_projected_points << std::endl;

				Projected_points boundary_clutter_projected;
				number_of_projected_points = m_ground_projector.project_with_planes(input, boundary_clutter, ground_plane, boundary_clutter_projected);
				log.out << "(7) Building's boundary clutter is projected. Number of projected points: " << number_of_projected_points << std::endl << std::endl;
				
				// Log points_exporter; points_exporter.export_projected_points_as_xyz("tmp/projected", building_boundaries_projected, m_default_path);

				
				// ----------------------------------

				// (7') Clean projected points by removing all points that lie far away from the center cluster of points.
				if (m_clean_projected_points) {

					std::cout << "(7') cleaning" << std::endl;
					m_preprocessor.set_scale(m_preprocessor_scale);

					auto number_of_removed_points = m_preprocessor.clean_projected_points(building_boundaries_projected, building_boundaries);
					log.out << "(7') Building's boundaries are cleaned. Number of removed points: " << number_of_removed_points << std::endl;

					number_of_removed_points = m_preprocessor.clean_projected_points(boundary_clutter_projected, boundary_clutter);
					log.out << "(7') Building's boundary clutter is cleaned. Number of removed points: " << number_of_removed_points << std::endl << std::endl;
				}


				// ----------------------------------

				// (8) Fit lines to the projected points in 2D.
				std::cout << "(8) line fitting" << std::endl;

				Lines lines;
				const auto number_of_fitted_lines = m_utils.fit_lines_to_projected_points(building_boundaries_projected, building_boundaries, lines);

				log.out << "(8) Lines are fitted. Number of fitted lines: " << number_of_fitted_lines << std::endl << std::endl;


				// ----------------------------------

				// (9) Find segments from the given lines.
				std::cout << "(9) creating segments" << std::endl;

				Segments segments;
				const auto number_of_segments = m_utils.create_segments_from_lines(building_boundaries_projected, building_boundaries, lines, segments);

				log.out << "(9) Segments are created. Number of created segments: " << number_of_segments << std::endl << std::endl;

				Log segments_exporter; segments_exporter.export_segments_as_obj("tmp/segments", segments, m_default_path);


				// ----------------------------------

				// (10) Apply 2D structuring algorithm.
				std::cout << "(10) 2d structuring" << std::endl;

				m_structuring = std::make_unique<Structuring_2>(building_boundaries_projected, building_boundaries, lines);
				
				m_structuring->set_epsilon(m_structuring_epsilon);
				m_structuring->save_log(m_structuring_log);
				m_structuring->resample(m_structuring_resample);

				const auto number_of_structured_segments = m_structuring->structure_point_set();

				log.out << "(10) 2D Structuring is applied. Number of structured segments: " << number_of_structured_segments << std::endl << std::endl;


				// ----------------------------------

				// (11) Compute constrained Delaunay triangulation of the structured points.
				std::cout << "(11) creating cdt" << std::endl;

				const Structured_points &structured_points = m_structuring_get_all_points ? m_structuring->get_structured_points() : m_structuring->get_segment_end_points();
				const Structured_labels &structured_labels = m_structuring_get_all_points ? m_structuring->get_structured_labels() : m_structuring->get_segment_end_labels();

				CDT cdt; 
				const auto number_of_faces = m_utils.compute_cdt(structured_points, structured_labels, cdt, 
													 m_add_cdt_clutter, boundary_clutter, boundary_clutter_projected, 
													 m_add_cdt_bbox, input);

				log.out << "(11) Constrained Delaunay triangulation of the structured points is built. Number of faces: " << number_of_faces << std::endl << std::endl;

				assert(!m_add_cdt_bbox); // visibility and graph cut do not work if bbox vertices are added to CDT!


				// ----------------------------------

				// (12) Convert 3D input to 2D input.			
				std::cout << "(12) converting 3d input into 2d input and setting face to points map" << std::endl;

				Container_2D input_2d; Face_points_map fp_map;
				const auto number_of_converted_points = m_utils.get_2d_input_and_face_points_map(cdt, input, input_2d, fp_map);

				log.out << "(12) 3D input is converted into 2D input and face to points map is set. Number of converted points: " << number_of_converted_points << std::endl << std::endl;


				// ----------------------------------

				// (13) Compute visibility (0 - outside or 1 - inside) for each triangle in CDT above.
				std::cout << "(13) visibility computation" << std::endl;

				m_visibility.save_info(m_visibility_save_info);
				m_visibility.set_approach(m_visibility_approach);
				m_visibility.set_method(m_visibility_method);
				m_visibility.set_number_of_samples(m_visibility_num_samples);

				const auto number_of_traversed_faces = m_visibility.compute(input_2d, cdt);
				log.out << "(13) Visibility is computed. Number of traversed faces: " << number_of_traversed_faces << std::endl << std::endl;

				// Log eps_saver_wp; eps_saver_wp.save_visibility_eps(cdt, input, structured_points); // works only with basic test
				
				Log eps_saver; eps_saver.save_visibility_eps(cdt);
				Log ply_vis_saver; ply_vis_saver.save_cdt_ply(cdt, "tmp/visibility", "in");


				// ----------------------------------

				// (14) Apply graph cut.
				std::cout << "(14) applying graph cut" << std::endl;

				m_graph_cut.save_info(m_graph_cut_save_info);
				m_graph_cut.set_alpha_parameter(m_graph_cut_alpha);
				m_graph_cut.set_beta_parameter(m_graph_cut_beta);
				m_graph_cut.set_gamma_parameter(m_graph_cut_gamma);

				m_graph_cut.max_flow(cdt);

				log.out << "(14) Graph cut is applied." << std::endl << std::endl;
				Log ply_cdt_in; ply_cdt_in.save_cdt_ply(cdt, "tmp/after_cut", "in");


				// ----------------------------------				

				// (15) From now on we handle each building separately.


				// (a) Split all buildings.
				std::cout << "(15 a) splitting buildings" << std::endl;

				Buildings buildings;
				const auto number_of_buildings = m_building_splitter.split(cdt, buildings);

				log.out << "(15 a) All buildings are found. Number of buildings: " << number_of_buildings << std::endl;


				// ----------------------------------

				if (m_use_boundaries) {

					// (b) Find building's walls.
					std::cout << "(15 b) finding boundaries" << std::endl;

					m_building_outliner.save_info(m_building_boundaries_save_internal_info);
					m_building_outliner.set_max_inner_iterations(m_building_boundaries_max_inner_iters);
					m_building_outliner.set_max_outer_iterations(m_building_boundaries_max_outer_iters);
					
					m_building_outliner.find_boundaries(cdt, buildings);

					log.out << "(15 b) All boundaries are found." << std::endl;
					Log log_bounds; log_bounds.save_buildings_info(cdt, buildings, "tmp/buildings_info_with_boundaries");
					
					// log_bounds.clear(); log_bounds.save_cdt_ply(cdt, "tmp/chosen_vertices"); // debugging info
				}


				// ----------------------------------

				// (c) Fit roof height for each building.
				std::cout << "(15 c) fitting roofs" << std::endl;
				switch (m_roof_fitter_type) {
					
					case Roof_fitter_type::MIN: 
						m_building_min_roof_fitter.fit_roof_heights(cdt, input, fp_map, ground_fitted_plane, buildings);
						break;

					case Roof_fitter_type::AVG: 
						m_building_avg_roof_fitter.fit_roof_heights(cdt, input, fp_map, ground_fitted_plane, buildings);
						break;

					case Roof_fitter_type::MAX: 
						m_building_max_roof_fitter.fit_roof_heights(cdt, input, fp_map, ground_fitted_plane, buildings);
						break;			

					default:
						assert(!"Wrong roof fitter type!");
						break;
				}
				
				log.out << "(15 c) All roofs are fitted." << std::endl << std::endl;
				Log log_roofs; log_roofs.save_buildings_info(cdt, buildings, "tmp/buildings_info_final");


				// ----------------------------------

				// (16) LOD0 reconstruction.

				m_lods.use_boundaries(m_use_boundaries);

				Ground ground_bbox;
				m_utils. template compute_ground_bbox<Ground, Ground_point>(input, ground_bbox);

				assert(!ground_bbox.empty());
				std::cout << "(16) reconstructing lod0" << std::endl;

				Mesh mesh_0; Mesh_facet_colors mesh_facet_colors_0;
				m_lods.reconstruct_lod0(cdt, buildings, ground_bbox, mesh_0, mesh_facet_colors_0);

				log.out << "(16) Final LOD0 is reconstructed." << std::endl << std::endl;
				Log lod_0_saver; lod_0_saver.save_mesh_as_ply(mesh_0, mesh_facet_colors_0, "LOD0");


				// ----------------------------------	
				
				// (17) LOD1 reconstruction.
				
				std::cout << "(17) reconstructing lod1" << std::endl;

				Mesh mesh_1; Mesh_facet_colors mesh_facet_colors_1;
				m_lods.reconstruct_lod1(cdt, buildings, ground_bbox, mesh_1, mesh_facet_colors_1);

				log.out << "(17) Final LOD1 is reconstructed." << std::endl;
				Log lod_1_saver; lod_1_saver.save_mesh_as_ply(mesh_1, mesh_facet_colors_1, "LOD1");


				// ----------------------------------

				// (END) Save log.
				std::cout << "... finishing\n" << std::endl;

				log.out << "\n\nFINISH EXECUTION";
				log.save("create_lods_ver_a");
			}
			
			// Version b.
			void create_lods_b() {

				// (START) Create log.
				std::cout << "\nstarting ..." << std::endl;
				Log log; log.out << "START EXECUTION\n\n\n";

				set_global_parameters();
				assert_global_parameters();


				// ----------------------------------

				// (1) Read data.
				std::cout << "(1) loading" << std::endl;

				Container_3D input;
				m_loader.get_data(m_default_path + ".ply", input);

				log.out << "(1) Data are loaded. Number of points: " << input.number_of_points() << std::endl << std::endl;


				// ----------------------------------

				// (2) Split data with respect to 4 different semantic labels.
				std::cout << "(2) selection" << std::endl;
				Indices ground_idxs, building_boundary_idxs;

				m_ground_selector.select_elements(input, std::back_inserter(ground_idxs));
				m_building_boundary_selector.select_elements(input, std::back_inserter(building_boundary_idxs));

				log.out << "(2 a) Ground is found. Number of elements: " 			   << ground_idxs.size() << std::endl;
				log.out << "(2 b) Building boundaries are found. Number of elements: " << building_boundary_idxs.size() << std::endl << std::endl;


				// ----------------------------------

				// (3) Create plane from the ground points.
				std::cout << "(3) ground plane fitting" << std::endl;

				const Plane_3 base_ground_plane = Plane_3(FT(0), FT(0), FT(1), FT(0));

				Plane_3 fitted_ground_plane;
				m_utils.fit_ground_plane(input, ground_idxs, fitted_ground_plane);
				
				log.out << "(3 a) Base ground plane is: " << base_ground_plane << std::endl;
				log.out << "(3 b) Data-fitted ground plane is: " << fitted_ground_plane << std::endl << std::endl;


				// ----------------------------------

				// (4) Project all boundary points onto the ground plane.
				std::cout << "(4) projecting" << std::endl;

				Projected_points building_boundaries_projected; 
				auto number_of_projected_points = m_ground_projector.project_with_indices(input, building_boundary_idxs, base_ground_plane, building_boundaries_projected);
				
				Log proj_saver; proj_saver.export_projected_points_as_xyz("tmp/projected_boundaries", building_boundaries_projected, m_default_path);
				log.out << "(4) Building's boundary points are projected. Number of projected points: " << number_of_projected_points << std::endl << std::endl;


				// ----------------------------------

				// (5) Compute constrained Delaunay triangulation of the input boundary points.
				std::cout << "(5) creating cdt" << std::endl;

				CDT cdt; 
				const auto number_of_faces = m_utils.compute_delaunay(building_boundaries_projected, cdt);

				log.out << "(5) Delaunay triangulation of the input boundary points is built. Number of faces: " << number_of_faces << std::endl << std::endl;


				// ----------------------------------

				// (6) Convert 3D input to 2D input.			
				std::cout << "(6) converting 3d input into 2d input and setting face to points map" << std::endl;

				Container_2D input_2d; Face_points_map fp_map;
				const auto number_of_converted_points = m_utils.get_2d_input_and_face_points_map(cdt, input, input_2d, fp_map);

				log.out << "(6) 3D input is converted into 2D input and face to points map is set. Number of converted points: " << number_of_converted_points << std::endl << std::endl;


				// ----------------------------------

				// (7) Compute visibility (0 - outside or 1 - inside) for each triangle in CDT above.
				std::cout << "(7) visibility computation" << std::endl;

				m_visibility.save_info(m_visibility_save_info);
				m_visibility.set_approach(m_visibility_approach);
				m_visibility.set_method(m_visibility_method);
				m_visibility.set_number_of_samples(m_visibility_num_samples);

				const auto number_of_traversed_faces = m_visibility.compute(input_2d, cdt);
				log.out << "(7) Visibility is computed. Number of traversed faces: " << number_of_traversed_faces << std::endl << std::endl;
				
				Log eps_saver; eps_saver.save_visibility_eps(cdt);
				Log ply_vis_saver; ply_vis_saver.save_cdt_ply(cdt, "tmp/visibility", "in");


				// ----------------------------------

				// (8) Apply graph cut.
				std::cout << "(8) applying graph cut" << std::endl;

				m_graph_cut.save_info(m_graph_cut_save_info);
				m_graph_cut.set_alpha_parameter(m_graph_cut_alpha);
				m_graph_cut.set_beta_parameter(m_graph_cut_beta);
				m_graph_cut.set_gamma_parameter(m_graph_cut_gamma);

				m_graph_cut.max_flow(cdt);

				log.out << "(8) Graph cut is applied." << std::endl << std::endl;
				Log ply_cdt_in; ply_cdt_in.save_cdt_ply(cdt, "tmp/after_cut", "in");


				// ----------------------------------				

				// (9) From now on we handle each building separately.


				// (a) Split all buildings.
				std::cout << "(9 a) splitting buildings" << std::endl;

				Buildings buildings;
				const auto number_of_buildings = m_building_splitter.split(cdt, buildings);

				log.out << "(9 a) All buildings are found. Number of buildings: " << number_of_buildings << std::endl;


				// ----------------------------------

				if (m_use_boundaries) {

					// (b) Find building's walls.
					std::cout << "(9 b) finding boundaries" << std::endl;

					m_building_outliner.save_info(m_building_boundaries_save_internal_info);
					m_building_outliner.set_max_inner_iterations(m_building_boundaries_max_inner_iters);
					m_building_outliner.set_max_outer_iterations(m_building_boundaries_max_outer_iters);
					
					m_building_outliner.find_boundaries(cdt, buildings);

					log.out << "(9 b) All boundaries are found." << std::endl;
					Log log_bounds; log_bounds.save_buildings_info(cdt, buildings, "tmp/buildings_info_with_boundaries");
				}


				// ----------------------------------

				// (c) Fit roof height for each building.
				std::cout << "(9 c) fitting roofs" << std::endl;
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
				
				log.out << "(9 c) All roofs are fitted." << std::endl << std::endl;
				Log log_roofs; log_roofs.save_buildings_info(cdt, buildings, "tmp/buildings_info_final");


				// ----------------------------------

				// (10 a) LOD0 reconstruction.

				m_lods.use_boundaries(m_use_boundaries);

				Ground ground_bbox;
				m_utils. template compute_ground_bbox<Ground, Ground_point>(input, ground_bbox);

				assert(!ground_bbox.empty());
				std::cout << "(10 a) reconstructing lod0" << std::endl;

				Mesh mesh_0; Mesh_facet_colors mesh_facet_colors_0;
				m_lods.reconstruct_lod0(cdt, buildings, ground_bbox, mesh_0, mesh_facet_colors_0);

				log.out << "(10 a) Final LOD0 is reconstructed." << std::endl;
				Log lod_0_saver; lod_0_saver.save_mesh_as_ply(mesh_0, mesh_facet_colors_0, "LOD0");


				// ----------------------------------	
				
				// (10 b) LOD1 reconstruction.
				
				std::cout << "(10 b) reconstructing lod1" << std::endl;

				Mesh mesh_1; Mesh_facet_colors mesh_facet_colors_1;
				m_lods.reconstruct_lod1(cdt, buildings, ground_bbox, mesh_1, mesh_facet_colors_1);

				log.out << "(10 b) Final LOD1 is reconstructed." << std::endl;
				Log lod_1_saver; lod_1_saver.save_mesh_as_ply(mesh_1, mesh_facet_colors_1, "LOD1");


				// ----------------------------------

				// (END) Save log.
				std::cout << "... finishing\n" << std::endl;

				log.out << "\n\nFINISH EXECUTION";
				log.save("create_lods_ver_b");
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
			
			Graph_cut m_graph_cut;
			Lods m_lods;

			Building_splitter m_building_splitter;
			Building_outliner m_building_outliner;
			
			Building_min_roof_fitter m_building_min_roof_fitter;
			Building_avg_roof_fitter m_building_avg_roof_fitter;
			Building_max_roof_fitter m_building_max_roof_fitter;

			std::unique_ptr<Structuring_2> m_structuring;


			// Global parameters.
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
			bool m_reject_planes;
			bool m_use_boundaries;

			
			// Assert default values of all global parameters.
			void assert_global_parameters() {
			
				assert(m_default_path != "default");

				assert(m_preprocessor_scale  != -FT(1));
				assert(m_structuring_epsilon != -FT(1));

				assert(!m_add_cdt_bbox);
				assert(m_visibility_num_samples != 0);

				assert(!(m_visibility_approach == Visibility_approach::FACE_BASED && m_visibility_method == Visibility_method::POINT_BASED_CLASSIFICATION));
				assert(!(m_visibility_approach == Visibility_approach::POINT_BASED && m_visibility_method == Visibility_method::FACE_BASED_COUNT));
				assert(!(m_visibility_approach == Visibility_approach::POINT_BASED && m_visibility_method == Visibility_method::FACE_BASED_NATURAL_NEIGHBOURS));
				assert(!(m_visibility_method == Visibility_method::FACE_BASED_BARYCENTRIC));

				assert(m_graph_cut_alpha != -FT(1));
				assert(m_graph_cut_beta  != -FT(1));
				assert(m_graph_cut_gamma != -FT(1));

				assert(m_building_boundaries_max_inner_iters != 0);
				assert(m_building_boundaries_max_outer_iters != 0);

				assert(m_max_reg_angle != -FT(1));
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


				// More important.
				m_reject_planes 			 = true;
				m_structuring_resample 	 	 = true;
				m_structuring_get_all_points = false;
				m_clean_projected_points 	 = false;
				
				m_visibility_approach = Visibility_approach::POINT_BASED;
				m_visibility_method   = Visibility_method::POINT_BASED_CLASSIFICATION;
				m_roof_fitter_type 	  = Roof_fitter_type::AVG;


				// The most important!
				const Main_test_data_type test_data_type = Main_test_data_type::BASIC;
				switch (test_data_type) {

					case Main_test_data_type::BASIC:
						set_basic_parameters();
						break;

					case Main_test_data_type::COMPLEX:
						set_complex_parameters();
						break;

					case Main_test_data_type::PARIS:
						set_paris_parameters();
						break;

					case Main_test_data_type::P10:
						set_p10_parameters();
						break;

					default:
						assert(!"Wrong test data!");
						break;
				}
			}

			void set_basic_parameters() {

				// To load basic parameters from stub, add stub to the loader class in LOD_traits!
				// Stub works with these basic parameters. Actually it is the same data set.

				// All main parameters are set below.
				m_default_path        	 = "/Users/danisimo/Documents/pipeline/data/basic_test/data";
				m_max_reg_angle          = 10.0;
				m_preprocessor_scale  	 = 2.0;
				m_structuring_epsilon 	 = 0.025; // the most important parameter!!!
				m_add_cdt_clutter     	 = true;
				m_visibility_num_samples = 200;
				m_graph_cut_alpha 		 = 1.0;
				m_graph_cut_beta 		 = 100000.0;
				m_graph_cut_gamma 		 = 1000.0;
				m_use_boundaries		 = true;
			}

			void set_complex_parameters() {

				// All main parameters are set below.
				m_default_path        	 = "/Users/danisimo/Documents/pipeline/data/complex_test/data_region_growing";
				m_max_reg_angle          = 10.0;
				m_preprocessor_scale  	 = 2.0;
				m_structuring_epsilon 	 = 0.0005; // the most important parameter!!!
				m_add_cdt_clutter     	 = false;
				m_visibility_num_samples = 200;
				m_graph_cut_alpha 		 = 1.0;
				m_graph_cut_beta 		 = 100000.0;
				m_graph_cut_gamma 		 = 1000.0;
				m_use_boundaries		 = false;
			}

			void set_paris_parameters() {

				// All main parameters are set below.
				m_default_path        	 = "/Users/danisimo/Documents/pipeline/data/paris_test/data_region_growing";
				m_max_reg_angle          = 10.0;
				m_preprocessor_scale  	 = 2.0;
				m_structuring_epsilon 	 = 0.1; // the most important parameter!!!
				m_add_cdt_clutter     	 = false;
				m_visibility_num_samples = 200;
				m_graph_cut_alpha 		 = 1.0;
				m_graph_cut_beta 		 = 100000.0;
				m_graph_cut_gamma 		 = 1000.0;
				m_use_boundaries		 = false;
			}

			void set_p10_parameters() {

				// All main parameters are set below.
				m_default_path        	 = "/Users/danisimo/Documents/pipeline/data/p10_test/data_region_growing";
				m_max_reg_angle          = 15.0;
				m_preprocessor_scale  	 = 2.0;
				m_structuring_epsilon 	 = 0.1; // the most important parameter!!!
				m_add_cdt_clutter     	 = false;
				m_visibility_num_samples = 200;
				m_graph_cut_alpha 		 = 1.0;
				m_graph_cut_beta 		 = 100000.0;
				m_graph_cut_gamma 		 = 1000.0;
				m_use_boundaries		 = false;
			}
		};
	}
}

#endif // CGAL_LEVEL_OF_DETAIL_BASE_H	