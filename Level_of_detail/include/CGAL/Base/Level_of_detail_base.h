#ifndef CGAL_LEVEL_OF_DETAIL_BASE_H
#define CGAL_LEVEL_OF_DETAIL_BASE_H

// STL includes.
#include <map>
#include <memory>
#include <string>
#include <iostream>
#include <vector>

// CGAL includes.
#include <CGAL/linear_least_squares_fitting_2.h>
#include <CGAL/linear_least_squares_fitting_3.h>
#include <CGAL/Triangulation_conformer_2.h>

// New CGAL includes.
#include <CGAL/Level_of_detail_enum.h>
#include <CGAL/Mylog/Mylog.h>

// Boost includes.
#include <boost/tuple/tuple.hpp>

namespace CGAL {

	namespace LOD {

		// Facade class that accumulates all necessary objects and operations
		// related to the level of detail (LOD) reconstruction.
		template<class LodTraits>
		class Level_of_detail_base {

		public:
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

			typedef typename Traits::Structuring_2 Structuring_2;
			typedef typename Traits::Visibility_2  Visibility_2;
			
			typedef typename Traits::CDT        		CDT;
			typedef typename CDT::Vertex_handle 		Vertex_handle;
			typedef typename CDT::Face_handle   		Face_handle;
			typedef typename CDT::Finite_edges_iterator Edge_iterator;
			typedef typename CDT::Finite_faces_iterator Face_iterator;

			typedef typename Traits::Lod_0 Lod_0;
			typedef typename Traits::Lod_1 Lod_1;

			typedef typename Traits::Mesh 			   Mesh;
			typedef typename Traits::Mesh_facet_colors Mesh_facet_colors;

			typedef typename Traits::Buildings Buildings;

			typedef typename Lod_1::Point  Ground_point;
			typedef typename Lod_1::Ground Ground;

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

			const std::string default_path = "/Users/danisimo/Documents/pipeline/data/basic_test_v3/";
			// const std::string default_path = "/Users/danisimo/Documents/pipeline/data/complex_test/";

			Level_of_detail_base(Traits traits = Traits()) : m_traits(traits) { } // Do I need to create an instance of these traits here?

			template<class OutputIterator>
			void create_lods(const std::string &, OutputIterator &&) {

				// (START) Create log.
				std::cout << "\nstarting ..." << std::endl;
				Log log; log.out << "START EXECUTION\n\n\n";


				// ----------------------------------

				// (1) Read data.
				std::cout << "(1) loading" << std::endl;

				Container_3D input;
				m_loader.get_data(default_path + "data.ply", input);

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

				Plane_3 ground_plane = Plane_3(0.0, 0.0, 1.0, 0.0); // use XY plane instead

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

				const auto number_of_regularized_planes = m_vertical_regularizer.regularize(building_boundaries, input, ground_plane);

				log.out << "(6) Building's nearly vertical planes are regularized. Number of regularized planes: " << number_of_regularized_planes <<
				", number of rejected planes: " << number_of_boundaries - building_boundaries.size() << std::endl << std::endl;

				// Log ply_saver; ply_saver.save_ply<Kernel, Container_3D>(input, "regularized", true);


				// ----------------------------------

				// (7) Project all vertical building's boundaries onto the ground plane.
				std::cout << "(7) projecting" << std::endl;

				Projected_points building_boundaries_projected; 
				auto number_of_projected_points = m_ground_projector.project(input, building_boundaries, ground_plane, building_boundaries_projected);
				log.out << "(7) Building's boundary planar points are projected. Number of projected points: " << number_of_projected_points << std::endl;

				Projected_points boundary_clutter_projected;
				number_of_projected_points = m_ground_projector.project(input, boundary_clutter, ground_plane, boundary_clutter_projected);
				log.out << "(7) Building's boundary clutter is projected. Number of projected points: " << number_of_projected_points << std::endl << std::endl;
				
				// Log points_exporter; points_exporter.export_projected_points_as_xyz("tmp/projected", building_boundaries_projected, default_path);

				
				// ----------------------------------

				// (7') Clean projected points by removing all points that lie far away from the center cluster of points.
				std::cout << "(7') cleaning" << std::endl;
				m_preprocessor.set_scale(FT(2));

				auto number_of_removed_points = m_preprocessor.clean_projected_points(building_boundaries_projected, building_boundaries);
				log.out << "(7') Building's boundaries are cleaned. Number of removed points: " << number_of_removed_points << std::endl;

				number_of_removed_points = m_preprocessor.clean_projected_points(boundary_clutter_projected, boundary_clutter);
				log.out << "(7') Building's boundary clutter is cleaned. Number of removed points: " << number_of_removed_points << std::endl << std::endl;


				// ----------------------------------

				// (8) Fit lines to the projected points in 2D.
				std::cout << "(8) line fitting" << std::endl;

				Lines lines;
				const auto number_of_fitted_lines = fit_lines_to_projected_points(building_boundaries_projected, building_boundaries, lines);

				log.out << "(8) Lines are fitted. Number of fitted lines: " << number_of_fitted_lines << std::endl << std::endl;


				// ----------------------------------

				// (9) Find segments from the given lines.
				std::cout << "(9) creating segments" << std::endl;

				Segments segments;
				const auto number_of_segments = create_segments_from_lines(building_boundaries_projected, building_boundaries, lines, segments);

				log.out << "(9) Segments are created. Number of created segments: " << number_of_segments << std::endl << std::endl;

				// Log segments_exporter; segments_exporter.export_segments_as_obj("tmp/segments", segments, default_path);


				// ----------------------------------

				// (10) Apply 2D structuring algorithm.
				std::cout << "(10) 2d structuring" << std::endl;

				m_structuring = std::make_unique<Structuring_2>(building_boundaries_projected, building_boundaries, lines);
				
				// epsilon: 0.025 basic test; 0.0005 complex test;
				m_structuring->set_epsilon(0.025);
				m_structuring->save_log(false);
				m_structuring->resample(true);

				const auto number_of_structured_segments = m_structuring->structure_point_set();

				log.out << "(10) 2D Structuring is applied. Number of structured segments: " << number_of_structured_segments << std::endl << std::endl;


				// ----------------------------------

				// (11) Compute constrained Delaunay triangulation of the structured points.
				std::cout << "(11) creating cdt" << std::endl;

				const Structured_points &structured_points = m_structuring->get_segment_end_points();
				const Structured_labels &structured_labels = m_structuring->get_segment_end_labels();

				// const Structured_points &structured_points = m_structuring->get_structured_points();
				// const Structured_points &structured_labels = m_structuring->get_structured_labels();

				CDT cdt; 

				// add_clutter: "true" for basic test and "false" for complex test
				const bool add_clutter = true; 
				const bool add_bbox    = false;

				const auto number_of_faces = compute_cdt(structured_points, structured_labels, cdt, 
														 add_clutter, boundary_clutter, boundary_clutter_projected, 
														 add_bbox, input);

				log.out << "(11) Constrained Delaunay triangulation of the structured points is built. Number of faces: " << number_of_faces << std::endl << std::endl;

				assert(!add_bbox); // visibility and graph cut do not work if bbox vertices are added to CDT!


				// ----------------------------------

				// (11') Convert 3D input to 2D input.			
				std::cout << "(11') converting 3d input into 2d input and setting face to points map" << std::endl;

				Container_2D input_2d; Face_points_map fp_map;
				const auto number_of_converted_points = get_2d_input_and_face_points_map(cdt, input, input_2d, fp_map);

				log.out << "(11') 3D input is converted into 2D input and face to points map is set. Number of converted points: " << number_of_converted_points << std::endl << std::endl;


				// ----------------------------------

				// (12) Compute visibility (0 - outside or 1 - inside) for each triangle in CDT above.
				std::cout << "(12) visibility computation" << std::endl;

				m_visibility.save_info(false);
				m_visibility.set_approach(Visibility_approach::POINT_BASED);
				m_visibility.set_method(Visibility_method::POINT_BASED_CLASSIFICATION);
				m_visibility.set_number_of_samples(200);

				const auto number_of_traversed_faces = m_visibility.compute(input_2d, cdt);
				log.out << "(12) Visibility is computed. Number of traversed faces: " << number_of_traversed_faces << std::endl << std::endl;

				Log eps_saver; eps_saver.save_visibility_eps(cdt, input, structured_points); // basic test
				// Log eps_saver; eps_saver.save_visibility_eps(cdt); 						 // complex test

				Log ply_vis_saver; ply_vis_saver.save_cdt_ply(cdt, "tmp/visibility", "in");


				// ----------------------------------

				// (13) Apply graph cut.
				std::cout << "(13) applying graph cut" << std::endl;

				m_lod_0.save_info(false);
				m_lod_0.set_alpha_parameter(FT(1));
				m_lod_0.set_beta_parameter(FT(100000));
				m_lod_0.set_gamma_parameter(FT(1000));

				m_lod_0.max_flow(cdt);

				log.out << "(13) Graph cut is applied." << std::endl << std::endl;
				Log ply_cdt_in; ply_cdt_in.save_cdt_ply(cdt, "tmp/after_cut", "in");


				// ----------------------------------				

				// (14) From now on we handle each building separately.


				// (a) Split all buildings.
				std::cout << "(14 a) splitting buildings" << std::endl;

				Buildings buildings;
				const auto number_of_buildings = m_building_splitter.split(cdt, buildings);

				log.out << "(14 a) All buildings are found. Number of buildings: " << number_of_buildings << std::endl;


				// ----------------------------------

				// (b) Find building's walls.
				std::cout << "(14 b) finding boundaries" << std::endl;

				m_building_outliner.save_info(false);
				m_building_outliner.set_max_inner_iterations(1000);
				m_building_outliner.set_max_outer_iterations(1000000);
				
				m_building_outliner.find_boundaries(cdt, buildings);

				log.out << "(14 b) All boundaries are found." << std::endl;
				Log log_bounds; log_bounds.save_buildings_info(cdt, buildings, "tmp/buildings_info_with_boundaries");


				// ----------------------------------

				// (c) Fit roof height for each building.
				std::cout << "(14 c) fitting roofs" << std::endl;

				const Roof_fitter_type roof_fitter_type = Roof_fitter_type::MAX;

				switch (roof_fitter_type) {
					
					case Roof_fitter_type::MIN: 
						m_building_min_roof_fitter.fit_roof_heights(cdt, input, fp_map, buildings);
						break;

					case Roof_fitter_type::AVG: 
						m_building_avg_roof_fitter.fit_roof_heights(cdt, input, fp_map, buildings);
						break;

					case Roof_fitter_type::MAX: 
						m_building_max_roof_fitter.fit_roof_heights(cdt, input, fp_map, buildings);
						break;

					default:
						assert(!"Wrong roof fitter type!");
						break;
				}
				
				log.out << "(14 c) All roofs are fitted." << std::endl << std::endl;
				Log log_roofs; log_roofs.save_buildings_info(cdt, buildings, "tmp/buildings_info_final");


				// ----------------------------------

				// (15) LOD0 reconstruction.

				Ground ground_bbox;
				compute_ground_bbox(input, ground_bbox);

				assert(!ground_bbox.empty());
				std::cout << "(15) reconstructing lod0" << std::endl;

				Mesh mesh_0; Mesh_facet_colors mesh_facet_colors_0;
				m_lod_1.reconstruct_lod0(cdt, buildings, ground_bbox, mesh_0, mesh_facet_colors_0);

				log.out << "(15) Final LOD0 is reconstructed." << std::endl << std::endl;
				Log lod_0_saver; lod_0_saver.save_mesh_as_ply(mesh_0, mesh_facet_colors_0, "LOD0");


				// ----------------------------------	
				
				// (16) LOD1 reconstruction.
				
				std::cout << "(16) reconstructing lod1" << std::endl;

				Mesh mesh_1; Mesh_facet_colors mesh_facet_colors_1;
				m_lod_1.reconstruct_lod1(cdt, buildings, ground_bbox, mesh_1, mesh_facet_colors_1);

				log.out << "(16) Final LOD1 is reconstructed." << std::endl;
				Log lod_1_saver; lod_1_saver.save_mesh_as_ply(mesh_1, mesh_facet_colors_1, "LOD1");


				// ----------------------------------

				// (END) Save log.
				std::cout << "... finishing\n" << std::endl;

				log.out << "\n\nFINISH EXECUTION";
				log.save("create_lods");
			}
			
		private:
			Traits       m_traits;
			Loader       m_loader;
			Preprocessor m_preprocessor;
			
			Building_boundary_selector m_building_boundary_selector;
			Building_interior_selector m_building_interior_selector;
			Clutter_selector           m_clutter_selector;
			Ground_selector            m_ground_selector;

			Vertical_regularizer	   m_vertical_regularizer;
			Ground_projector 		   m_ground_projector;
			Visibility_2 			   m_visibility;
			
			Lod_0 m_lod_0;
			Lod_1 m_lod_1;

			Building_splitter m_building_splitter;
			Building_outliner m_building_outliner;
			
			Building_min_roof_fitter m_building_min_roof_fitter;
			Building_avg_roof_fitter m_building_avg_roof_fitter;
			Building_max_roof_fitter m_building_max_roof_fitter;

			std::unique_ptr<Structuring_2> m_structuring;
			
			// Not efficient since I need to copy all ground points.
			void fit_ground_plane(const Container_3D &input, const Indices &ground_idxs, Plane_3 &ground_plane) const {

				std::vector<Point_3> tmp_ground(ground_idxs.size());
				for (size_t i = 0; i < ground_idxs.size(); ++i) tmp_ground[i] = input.point(ground_idxs[i]);

				CGAL::linear_least_squares_fitting_3(tmp_ground.begin(), tmp_ground.end(), ground_plane, CGAL::Dimension_tag<0>()); 
			}

			// Not efficient since I need to copy all ground points.
			int fit_lines_to_projected_points(const Projected_points &points, const Planes &planes, Lines &lines) const {

				auto number_of_fitted_lines = 0;
				std::vector<Point_2> tmp_points;

				lines.clear();
				lines.resize(planes.size());

				for (Plane_iterator it = planes.begin(); it != planes.end(); ++it, ++number_of_fitted_lines) {
					const auto num_points = (*it).second.size();

					tmp_points.clear();
					tmp_points.resize(num_points);

					for (size_t i = 0; i < num_points; ++i) {
						const auto index = (*it).second[i];
						tmp_points[i] = points.at(index);
					}

					CGAL::linear_least_squares_fitting_2(tmp_points.begin(), tmp_points.end(), lines[number_of_fitted_lines], CGAL::Dimension_tag<0>());
				}
				return number_of_fitted_lines;
			}

			// It may have precision problems: e.g. Line = (1 1.0e+24) -- (1 1.0e+24) or nan!
			int create_segments_from_lines(const Projected_points &points, const Planes &planes, const Lines &lines, Segments &segments) const {

				segments.clear();
				segments.resize(lines.size());

				std::cout.precision(20);

				auto number_of_segments = 0;
				for (Plane_iterator it = planes.begin(); it != planes.end(); ++it, ++number_of_segments) {

					const auto num_points = (*it).second.size();
					const auto big = +100000.0;

					auto minx = big, miny = big;
					auto maxx = -minx, maxy = -miny;				

					for (size_t i = 0; i < num_points; ++i) {

						const auto index = (*it).second[i];
						const auto point = points.at(index);

						const auto projected = project(lines[number_of_segments], point);

						minx = CGAL::min(minx, projected.x());
						maxx = CGAL::max(maxx, projected.x());

						miny = CGAL::min(miny, projected.y());
						maxy = CGAL::max(maxy, projected.y());
					}
					segments[number_of_segments] = Segment_2(Point_2(minx, miny), Point_2(maxx, maxy));

					auto v1 = lines[number_of_segments].to_vector();
					auto v2 = segments[number_of_segments].to_vector();

					// Rotate segments if needed.
					const auto eps = 1.0e-6;
					if ((v1.y() < 0.0 && v2.y() >= 0.0 && std::fabs(v1.y() - v2.y()) > eps) ||
						(v2.y() < 0.0 && v1.y() >= 0.0 && std::fabs(v1.y() - v2.y()) > eps) ||
						(v1.x() < 0.0 && v2.x() >= 0.0 && std::fabs(v1.x() - v2.x()) > eps) ||
						(v2.x() < 0.0 && v1.x() >= 0.0 && std::fabs(v1.x() - v2.x()) > eps)) {

						segments[number_of_segments] = Segment_2(Point_2(minx, maxy), Point_2(maxx, miny));
					}
				}

				assert(static_cast<size_t>(number_of_segments) == lines.size());
				return number_of_segments;
			}

			// My custom function to handle precision problems when projecting points.
			Point_2 project(const Line_2 &line, const Point_2 &p) const {

				const auto a = line.point(0);
				const auto b = line.point(1);

				const auto projected = a + CGAL::scalar_product(p - a, b - a) / CGAL::scalar_product(b - a, b - a) * (b - a);
				
				if (std::isnan(projected.x()) || std::isnan(projected.y())) return line.projection(p);
				else return projected;
			}

			// BE CAREFUL: THIS THING CAN INSERT NEW POINTS!
			int compute_cdt(const Structured_points &points, const Structured_labels &labels, CDT &cdt, 
							const bool add_clutter, const Boundary_data &boundary_clutter, const Projected_points &boundary_clutter_projected,
							const bool add_bbox   , const Container_3D  &input) const {

				Log log;
				auto number_of_faces = -1;

				assert(!points.empty() && !labels.empty());
				assert(points.size() == labels.size());


				// Add all structured segments/points with the corresponding labels.
				std::vector<std::vector<Vertex_handle> > vhs(points.size());

				// Insert points with labels.
				for (size_t i = 0; i < points.size(); ++i) {
					assert(points[i].size() == labels[i].size());

					vhs[i].resize(points[i].size());
					for (size_t j = 0; j < points[i].size(); ++j) {

						vhs[i][j] = cdt.insert(points[i][j]);
						vhs[i][j]->info().label = labels[i][j];
					}
				}

				// Insert constraints.
				for (size_t i = 0; i < points.size(); ++i)
					for (size_t j = 0; j < points[i].size() - 1; ++j)
						cdt.insert_constraint(vhs[i][j], vhs[i][j + 1]);


				// Add clutter.
				if (add_clutter) {
					for (typename Boundary_data::const_iterator it = boundary_clutter.begin(); it != boundary_clutter.end(); ++it) {
						const size_t num_clutter_points = (*it).second.size();

						for (size_t i = 0; i < num_clutter_points; ++i) {
							const auto point_index = (*it).second[i];

							Vertex_handle vh = cdt.insert(boundary_clutter_projected.at(point_index));
							vh->info().label = Structured_label::CLUTTER;
						}
					}
				}


				// Add bounding box.
				if (add_bbox) {
					
					std::vector<Point_2> bbox;
					compute_bounding_box(input, bbox);

					std::vector<Vertex_handle> bhs(bbox.size());

					for (size_t i = 0; i < bbox.size(); ++i) bhs[i] = cdt.insert(bbox[i]);
					for (size_t i = 0; i < bbox.size(); ++i) {

						const size_t ip = (i + 1) % bbox.size();
						cdt.insert_constraint(bhs[i], bhs[ip]);
					}
				}


				// Create CDT.
				CGAL::make_conforming_Delaunay_2(cdt);
				number_of_faces = cdt.number_of_faces();


				// Save CDT.
				log.save_cdt_obj(cdt, "tmp/cdt");
				return number_of_faces;
			}

			void compute_ground_bbox(const Container_3D &input, Ground &ground_bbox) const {

				ground_bbox.clear();
				ground_bbox.resize(4);

				std::vector<Point_2> bbox;
				compute_bounding_box(input, bbox);

				assert(bbox.size() == 4);
				ground_bbox[0] = Ground_point(bbox[0].x(), bbox[0].y(), FT(0));
				ground_bbox[1] = Ground_point(bbox[1].x(), bbox[1].y(), FT(0));
				ground_bbox[2] = Ground_point(bbox[2].x(), bbox[2].y(), FT(0));
				ground_bbox[3] = Ground_point(bbox[3].x(), bbox[3].y(), FT(0));
			}

			void compute_bounding_box(const Container_3D &input, std::vector<Point_2> &bbox) const {

				bbox.clear();
				bbox.resize(4);

				const FT big_value = FT(1000000); // change it

				FT minx =  big_value, miny =  big_value;
				FT maxx = -big_value, maxy = -big_value;

				for (typename Container_3D::const_iterator it = input.begin(); it != input.end(); ++it) {
					const Point_3 &p = input.point(*it);

					const FT x = p.x();
					const FT y = p.y();

					minx = CGAL::min(minx, x);
					miny = CGAL::min(miny, y);

					maxx = CGAL::max(maxx, x);
					maxy = CGAL::max(maxy, y);
				}

				bbox[0] = Point_2(minx, miny);
				bbox[1] = Point_2(maxx, miny);
				bbox[2] = Point_2(maxx, maxy);
				bbox[3] = Point_2(minx, maxy);
			}

			int get_2d_input_and_face_points_map(const CDT &cdt, const Container_3D &input_3d, Container_2D &input_2d, Face_points_map &fp_map) const {

				Log log;

				typename CDT::Locate_type locate_type;
				int locate_index_stub = -1;

				Label_map class_labels;
				boost::tie(class_labels, boost::tuples::ignore) = input_3d.template property_map<Label>("label");

				input_2d.clear();
				input_2d.resize(input_3d.number_of_points());

				size_t point_index = 0;
				for (typename Container_3D::const_iterator it = input_3d.begin(); it != input_3d.end(); ++it, ++point_index) {

					const Point_index pi = *it;

					const Point_3 &p = input_3d.point(pi);
					const Point_2 &q = Point_2(p.x(), p.y());

					input_2d[point_index] = std::make_pair(q, class_labels[pi]);
					log.out << q << " " << 0 << std::endl;

					const Face_handle face_handle = cdt.locate(q, locate_type, locate_index_stub);
					if (locate_type == CDT::FACE || locate_type == CDT::EDGE || locate_type == CDT::VERTEX) fp_map[face_handle].push_back(pi);
				}

				log.save("tmp/input_2d", ".xyz");
				return point_index;
			}

			void add_bbox_to(const CDT &cdt, const Container_3D &input, CDT &cdt_with_bbox) {

				CDT tmp = cdt;

				for (Edge_iterator eit = tmp.finite_edges_begin(); eit != tmp.finite_edges_end(); ++eit)
					if (is_boundary_edge(tmp, eit)) insert_constraint(tmp, eit);

				cdt_with_bbox = tmp;

				std::vector<Point_2> bbox;
				compute_bounding_box(input, bbox);

				std::vector<Vertex_handle> bhs(bbox.size());

				for (size_t i = 0; i < bbox.size(); ++i) bhs[i] = cdt_with_bbox.insert(bbox[i]);
				for (size_t i = 0; i < bbox.size(); ++i) {

					const size_t ip = (i + 1) % bbox.size();
					cdt_with_bbox.insert_constraint(bhs[i], bhs[ip]);
				}
			}

			bool is_boundary_edge(const CDT &cdt, const Edge_iterator &edge_handle) const {

				const int vertex_index = edge_handle->second;

				const Face_handle face_1 = edge_handle->first;
				const Face_handle face_2 = face_1->neighbor(vertex_index);

				if (cdt.is_infinite(face_1) || cdt.is_infinite(face_2)) return true;
				return false;
			}

			void insert_constraint(CDT &cdt, const Edge_iterator &edge_handle) const {

				const Face_handle face = edge_handle->first;
				const int vertex_index = edge_handle->second;

				const Vertex_handle &source = face->vertex(cdt.ccw(vertex_index));
				const Vertex_handle &target = face->vertex( cdt.cw(vertex_index));

				cdt.insert_constraint(source, target);
			}
		};
	}
}

#endif // CGAL_LEVEL_OF_DETAIL_BASE_H	