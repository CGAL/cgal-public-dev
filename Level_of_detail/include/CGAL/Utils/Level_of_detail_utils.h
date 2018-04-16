#ifndef CGAL_LEVEL_OF_DETAIL_UTILS_H
#define CGAL_LEVEL_OF_DETAIL_UTILS_H

#if defined(WIN32) || defined(_WIN32) 
#define PSR "\\" 
#else 
#define PSR "/" 
#endif 

// STL includes.
#include <map>
#include <cassert>
#include <vector>
#include <memory>
#include <string>

// CGAL includes.
#include <CGAL/Triangulation_conformer_2.h>
#include <CGAL/Barycentric_coordinates_2/Segment_coordinates_2.h>
#include <CGAL/linear_least_squares_fitting_2.h>
#include <CGAL/linear_least_squares_fitting_3.h>
#include <CGAL/utils.h>
#include <CGAL/Random.h>
#include <CGAL/number_utils.h>
#include <CGAL/constructions_d.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Kd_tree.h>
#include <CGAL/Search_traits_2.h>
#include <CGAL/Search_traits_adapter.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Fuzzy_sphere.h>
#include <CGAL/property_map.h>

// New CGAL includes.
#include <CGAL/Mylog/Mylog.h>
#include <CGAL/Level_of_detail_enum.h>
#include <CGAL/Utils/Level_of_detail_utils_simple.h>
#include <CGAL/Utils/Level_of_detail_uniform_sample_generator.h>

// Boost includes.
#include <boost/tuple/tuple.hpp>

namespace BC = CGAL::Barycentric_coordinates;

namespace CGAL {

	namespace LOD {

		template<class KernelTraits, class InputContainer, class InputCDT>
		class Level_of_detail_utils {

		public:
			typedef KernelTraits   Kernel;
			typedef InputContainer Input;
			typedef InputCDT 	   CDT;

			typedef typename Kernel::FT 	   FT;
			typedef typename Kernel::Point_2   Point_2;
			typedef typename Kernel::Point_3   Point_3;
			typedef typename Kernel::Plane_3   Plane_3;
			typedef typename Kernel::Vector_3  Vector_3;
			typedef typename Kernel::Line_2    Line_2;
			typedef typename Kernel::Line_3    Line_3;
			typedef typename Kernel::Segment_2 Segment_2;

			typedef typename CDT::Vertex_handle 		   Vertex_handle;
			typedef typename CDT::Face_handle 			   Face_handle;
			typedef typename CDT::Finite_edges_iterator    Edge_iterator;
			typedef typename CDT::Finite_vertices_iterator Vertex_iterator;
			typedef typename CDT::Finite_faces_iterator    Face_iterator;

			typedef Level_of_detail_utils_simple<Kernel> Simple_utils;

			// Extra.
			using Label       = int;
			using Label_map   = typename Input:: template Property_map<Label>;
			using Point_index = typename Input::Index;

			using Segments = std::vector<Segment_2>;
			using Log 	   = CGAL::LOD::Mylog;

			typename Kernel::Compute_squared_distance_2 	  squared_distance;
			typename Kernel::Compute_squared_length_3 		  squared_length;
			typename Kernel::Compute_scalar_product_3 		  dot_product;
			typename Kernel::Construct_cross_product_vector_3 cross_product;


		private:
			Simple_utils m_simple_utils;

			void create_samples(const Point_2 &source, const Point_2 &target, std::vector<Point_2> &samples) const {
				samples.clear();

				const size_t num_steps = 5;

				assert(num_steps > 1);
				samples.resize(num_steps - 1);

				for (size_t i = 1; i < num_steps; ++i) {

					const FT t = static_cast<FT>(i) / static_cast<FT>(num_steps);
					const FT s = FT(1) - t;

					const FT x = s * source.x() + t * target.x();
					const FT y = s * source.y() + t * target.y();

					samples[i-1] = Point_2(x, y);
				}
			}

		public:

			///////////////////////////
			// Mesh functions!

			template<class Mesh>
			void translate_mesh_along_z(Mesh &mesh, const FT target_height) const {

				using Polyhedron_vertex_iterator = typename Mesh::Vertex_iterator;
				for (Polyhedron_vertex_iterator vit = mesh.vertices_begin(); vit != mesh.vertices_end(); ++vit) {

					Point_3 &point = vit->point();
					point = Point_3(point.x(), point.y(), point.z() + target_height);
				}
			}

			///////////////////////////
			// Triangulation functions!

			template<class Data_structure>
			void update_constraints(Data_structure &data_structure, const Segments &segments, const bool make_silent) const {
				
				using Containers  = typename Data_structure::Containers;		
				using Constraints = typename Data_structure::Container::Constraints;

				using Container = typename Data_structure::Container;
				using Polygon   = typename Container::Polygon;

				using Polygon_vertex_iterator = typename Polygon::Vertex_const_iterator;

				Containers &containers = data_structure.containers();
				assert(containers.size() > 0);

				Segments updated_constraints;
				std::vector<Point_2> vertices;

				for (size_t i = 0; i < containers.size(); ++i) {
					
					const Polygon &polygon   = containers[i].polygon;
					Constraints &constraints = containers[i].constraints;

					const size_t num_edges = std::distance(polygon.vertices_begin(), polygon.vertices_end());

                	constraints.clear();
                	constraints.resize(num_edges, false);

					const size_t num_vertices = num_edges;
					vertices.resize(num_vertices);
					
					size_t j = 0;
                	for (Polygon_vertex_iterator vit = polygon.vertices_begin(); vit != polygon.vertices_end(); ++vit, ++j) 
						vertices[j] = *vit;

					for (j = 0; j < vertices.size(); ++j) {
						const size_t jp = (j + 1) % vertices.size();

						const Point_2 &p1 = vertices[j];
						const Point_2 &p2 = vertices[jp];

						constraints[j] = get_new_constraint(p1, p2, segments);

						if (!make_silent && constraints[j]) 
							updated_constraints.push_back(Segment_2(p1, p2));
					}
				}

				if (!make_silent) {
					const std::string stub = "";

                    Log segments_exporter; 
					segments_exporter.export_segments_as_obj("tmp" + std::string(PSR) + "updated_constraints", updated_constraints, stub);
				}
			}

			bool get_new_constraint(const Point_2 &p1, const Point_2 &p2, const Segments &segments) const {
				
				assert(segments.size() > 0);
				for (size_t i = 0; i < segments.size(); ++i) {
					
					const Segment_2 &segment = segments[i];
					if (is_constrained(p1, p2, segment)) return true;
				}

				return false;
			}

			bool is_constrained(const Point_2 &p1, const Point_2 &p2, const Segment_2 &segment) const {

				Segment_2 edge(p1, p2);

				const Point_2 &source = segment.source();
				const Point_2 &target = segment.target();

				Line_2 line(source, target);
				
				const Point_2 pr1 = line.projection(p1);
				const Point_2 pr2 = line.projection(p2);

				const FT eps = FT(3);

				if (squared_distance(p1, pr1) > eps * eps) return false;
				if (squared_distance(p2, pr2) > eps * eps) return false;

				const FT tol = -FT(1) / FT(10);
				
				std::pair<FT, FT> bc = BC::compute_segment_coordinates_2(source, target, p1, Kernel());
				const bool state1 = bc.first > tol && bc.second > tol;

				bc = BC::compute_segment_coordinates_2(source, target, p2, Kernel());
				const bool state2 = bc.first > tol && bc.second > tol;

				bc = BC::compute_segment_coordinates_2(p1, p2, source, Kernel());
				const bool state3 = bc.first > tol && bc.second > tol;

				bc = BC::compute_segment_coordinates_2(p1, p2, target, Kernel());
				const bool state4 = bc.first > tol && bc.second > tol;

				if ( (state1 && state2) || (state3 && state4) ) return true;
				return false;
			}

			void compute_face_barycentre(const Face_iterator &fit, Point_2 &barycentre) const {

				FT x = FT(0), y = FT(0);
				for (size_t i = 0; i < 3; ++i) {

					x += fit->vertex(i)->point().x();
					y += fit->vertex(i)->point().y();
				}

				x /= FT(3);
				y /= FT(3);

				barycentre = Point_2(x, y);
			}

			template<class Data_structure>
			void update_labels(CDT &cdt, const Data_structure &data_structure) const {
				
				using Containers = typename Data_structure::Containers;
				using Container  = typename Data_structure::Container;
				using Polygon 	 = typename Container::Polygon;
				
				const Containers &containers = data_structure.containers();
				assert(containers.size() > 0);

				Point_2 barycentre;
				for (Face_iterator fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit) {
					compute_face_barycentre(fit, barycentre);

					for (size_t i = 0; i < containers.size(); ++i) {
						const Polygon &polygon = containers[i].polygon;

						if (polygon.has_on_bounded_side(barycentre)) {
							fit->info().in       = containers[i].inside;
							fit->info().in_color = containers[i].colour;

							break;
						}
					}
				}
			}

			void clean_labels(CDT &cdt) const {
				
				const size_t num_faces = cdt.number_of_faces();

				std::vector<FT> 		 labels(num_faces);
				std::vector<CGAL::Color> colors(num_faces);

				size_t i = 0;
				for (Face_iterator fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit, ++i) {
					const FT label = fit->info().in;

					const Face_handle fh1 = fit->neighbor(0);
					const Face_handle fh2 = fit->neighbor(1);
					const Face_handle fh3 = fit->neighbor(2);

					const FT label1 = fh1->info().in;
					const FT label2 = fh2->info().in;
					const FT label3 = fh3->info().in;

					// Wrong labels.
					if ((label == FT(1) && label1 == FT(0) && label2 == FT(0) && label3 == FT(0)) ||
						(label == FT(0) && label1 == FT(1) && label2 == FT(1) && label3 == FT(1)) ){
						
						labels[i] = label1;
						colors[i] = fh1->info().in_color;

						continue;
					}

					// This face is too thin.
					/*
					const Point_2 &v1 = fit->vertex(0)->point();
					const Point_2 &v2 = fit->vertex(1)->point();
					const Point_2 &v3 = fit->vertex(2)->point();

					const typename Kernel::Triangle_2 tri(v1, v2, v3);
					if (tri.area() < FT(1)) {

						labels[i] = FT(0);
						colors[i] = CGAL::Color(255, 55, 55);

						continue;
					} */

					// Other faces.
					labels[i] = label;
					colors[i] = fit->info().in_color;
				}

				i = 0;
				for (Face_iterator fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit, ++i) {

					fit->info().in 		 = labels[i];
					fit->info().in_color = colors[i];
				}
			}

			template<class Data_structure, class Projected_points>
			void compute_cdt(CDT &cdt, const Data_structure &data_structure, const Projected_points &boundary_clutter_projected, const bool add_clutter, const bool make_silent) const {

				using Containers = typename Data_structure::Containers;
				using Container  = typename Data_structure::Container;
				using Polygon    = typename Container::Polygon;

				using Polygon_vertex_iterator = typename Polygon::Vertex_const_iterator;
				const Containers &containers = data_structure.containers();
				
				assert(containers.size() > 0);
				std::vector< std::vector<Vertex_handle> > vhs(containers.size());

				// Insert points.
				for (size_t i = 0; i < containers.size(); ++i) {
					const Polygon &polygon = containers[i].polygon;

					vhs[i].resize(std::distance(polygon.vertices_begin(), polygon.vertices_end()));

					size_t j = 0;
					for (Polygon_vertex_iterator vit = polygon.vertices_begin(); vit != polygon.vertices_end(); ++vit, ++j) {
						
						vhs[i][j] = cdt.insert(*vit);
						vhs[i][j]->info().label = Structured_label::CLUTTER;
					}
				}

				// Add clutter.
				if (add_clutter) {

					using Point_iterator = typename Projected_points::const_iterator;

					for (Point_iterator pit = boundary_clutter_projected.begin(); pit != boundary_clutter_projected.end(); ++pit) {
						const auto point = (*pit).second;

						Vertex_handle vh = cdt.insert(point);
						vh->info().label = Structured_label::CLUTTER;
					}
				}

				// Insert constraints.
				for (size_t i = 0; i < vhs.size(); ++i){
					for (size_t j = 0; j < vhs[i].size(); ++j) {
						const size_t jp = (j + 1) % vhs[i].size();
						
						if (vhs[i][j] != vhs[i][jp])
							cdt.insert_constraint(vhs[i][j], vhs[i][jp]);
					}
				}

				// Update labels.
				update_labels(cdt, data_structure);

				// Clean labels.
				clean_labels(cdt);

				// Save CDT.
				Log exporter;
				if (!make_silent) exporter.save_cdt_ply(cdt, "tmp" + std::string(PSR) + "lod_0_1" + std::string(PSR) + "cdt");
			}

			int compute_cdt(CDT &cdt, const Segments &segments, const bool silent, const bool sample = false) const {

				Log log;
				auto number_of_faces = -1;

				assert(segments.size() > 0);
				std::vector<Vertex_handle> vhs(segments.size() * 2);

				
				// Insert points with labels.
				std::vector<Point_2> samples;

				size_t count = 0;
				for (size_t i = 0; i < segments.size(); ++i) {

					const Point_2 &source = segments[i].source();
					const Point_2 &target = segments[i].target();

					vhs[count] = cdt.insert(source);
					vhs[count]->info().label = Structured_label::CLUTTER;
					++count;

					vhs[count] = cdt.insert(target);
					vhs[count]->info().label = Structured_label::CLUTTER;
					++count;

					if (sample) {
						create_samples(source, target, samples);
						for (size_t j = 0; j < samples.size(); ++j) {

							Vertex_handle vh = cdt.insert(samples[j]);
							vh->info().label = Structured_label::CLUTTER;
						}
					}
				}


				// Insert constraints.
				for (size_t i = 0; i < vhs.size(); i += 2)
					if (vhs[i] != vhs[i+1])
						cdt.insert_constraint(vhs[i], vhs[i+1]);

				
				if (!silent) {
					Log logex;
					logex.export_segments_as_obj("tmp" + std::string(PSR) + "cdt_constraints", segments, "stub");
				}


				// Get number of faces.
				number_of_faces = cdt.number_of_faces();


				// Save CDT.
				if (!silent) log.save_cdt_obj(cdt, "tmp" + std::string(PSR) + "cdt");
				return number_of_faces;
			}

			// BE CAREFUL: THIS THING CAN INSERT NEW POINTS!
			template<class Structured_points, class Structured_labels, class Structured_anchors, class Boundary_data, class Projected_points>
			int compute_cdt(const Structured_points &points, const Structured_labels &labels, const Structured_anchors &anchors, const FT adjacency_value, CDT &cdt, 
							const bool add_clutter, const Boundary_data &, const Projected_points &boundary_clutter_projected, const bool silent) const {

				Log log;
				auto number_of_faces = -1;

				assert(points.size() == labels.size());
				cdt.clear();

				assert(points.size() + boundary_clutter_projected.size() > 2);

				// Add all structured segments/points with the corresponding labels.
				std::vector<std::vector<Vertex_handle> > vhs(points.size());


				// Insert points with labels.
				Point_2 tmp(FT(1000000000), FT(1000000000));
				for (size_t i = 0; i < points.size(); ++i) {
					if (points[i].size() < 2) continue;

					assert(points[i].size() == labels[i].size());
					vhs[i].resize(points[i].size());

					for (size_t j = 0; j < points[i].size(); ++j) {
						if (tmp == points[i][j]) continue;

						vhs[i][j] = cdt.insert(points[i][j]);
						vhs[i][j]->info().label = labels[i][j];

						assert(labels[i][j] != Structured_label::CLUTTER);
						tmp = points[i][j];
					}
				}


				// Insert constraints.
				std::vector<Segment_2> constraints;
				for (size_t i = 0; i < points.size(); ++i) {
					if (points[i].size() < 2) continue;

					for (size_t j = 0; j < points[i].size() - 1; ++j) {	
						if (vhs[i][j] != Vertex_handle() && vhs[i][j + 1] != Vertex_handle() && is_valid_segment(points[i][j], points[i][j + 1])) {
							
							if (!silent) constraints.push_back(Segment_2(points[i][j], points[i][j + 1]));
							cdt.insert_constraint(vhs[i][j], vhs[i][j + 1]);
						}
					}
				}

				assert(adjacency_value > FT(0));
				for (size_t i = 0; i < points.size(); ++i) {
					if (points[i].size() < 2) continue;
					
					for (size_t j = 0; j < points[i].size(); ++j) {

						if (vhs[i][j] != Vertex_handle() && vhs[i][j]->info().label == Structured_label::CORNER) {
							assert(anchors[i][j].size() == 2);

							const int ind_a = anchors[i][j][0];
							const int ind_b = anchors[i][j][1];

							assert(ind_a >= 0 && ind_b >= 0);

							const int a_k = find_closest_point(points[i][j], points, ind_a);
							const int b_k = find_closest_point(points[i][j], points, ind_b);

							const FT eps = FT(1) / FT(100000);
							assert(a_k >= 0 && b_k >= 0);

							const FT dist_a = static_cast<FT>(CGAL::sqrt(CGAL::to_double(squared_distance(points[i][j], points[ind_a][a_k]))));
							const FT scale  = FT(1);

							if (dist_a > eps && dist_a < scale * adjacency_value) {
								if (!silent) constraints.push_back(Segment_2(points[i][j], points[ind_a][a_k]));

								assert(static_cast<int>(i) != ind_a);
								cdt.insert_constraint(vhs[i][j], vhs[ind_a][a_k]);
							}

							const FT dist_b = static_cast<FT>(CGAL::sqrt(CGAL::to_double(squared_distance(points[i][j], points[ind_b][b_k]))));

							if (dist_b > eps && dist_b < scale * adjacency_value) {
								if (!silent) constraints.push_back(Segment_2(points[i][j], points[ind_b][b_k]));

								assert(static_cast<int>(i) != ind_b);
								cdt.insert_constraint(vhs[i][j], vhs[ind_b][b_k]);
							}
						}
					}
				}

				if (!silent) {
					Log logex;
					logex.export_segments_as_obj("tmp" + std::string(PSR) + "cdt_constraints", constraints, "stub");
				}


				// Correct all wrong labels.
				for (Vertex_iterator vit = cdt.finite_vertices_begin(); vit != cdt.finite_vertices_end(); ++vit)
					if(vit->info().label == Structured_label::CLUTTER)
						vit->info().label = Structured_label::LINEAR;


				// Add clutter.
				if (add_clutter) {
					using Point_iterator = typename Projected_points::const_iterator;

					for (Point_iterator pit = boundary_clutter_projected.begin(); pit != boundary_clutter_projected.end(); ++pit) {
						const auto point = (*pit).second;

						Vertex_handle vh = cdt.insert(point);
						vh->info().label = Structured_label::CLUTTER;
					}

				} else {
					for (Vertex_iterator vit = cdt.finite_vertices_begin(); vit != cdt.finite_vertices_end(); ++vit)
						assert(vit->info().label != Structured_label::CLUTTER);
				}


				// Create CDT.
				// CGAL::make_conforming_Delaunay_2(cdt);
				number_of_faces = cdt.number_of_faces();


				// Correct all wrong labels again.
				if (!add_clutter) {
					for (Vertex_iterator vit = cdt.finite_vertices_begin(); vit != cdt.finite_vertices_end(); ++vit)
						if(vit->info().label == Structured_label::CLUTTER)
							vit->info().label = Structured_label::LINEAR;
				}


				// Save CDT.
				if (!silent) log.save_cdt_obj(cdt, "tmp" + std::string(PSR) + "cdt");
				return number_of_faces;
			}

			bool is_valid_segment(const Point_2 &a, const Point_2 &b) const {

				const FT eps = FT(1) / FT(100000);
				if (static_cast<FT>(CGAL::sqrt(CGAL::to_double(squared_distance(a, b)))) < eps) return false;

				return true;
			}

			template<class Structured_points>
			int find_closest_point(const Point_2 &corner, const Structured_points &points, const int segment_index) const {

				std::vector<FT> dist(points[segment_index].size());
				for (size_t i = 0; i < points[segment_index].size(); ++i) {

					dist[i] = static_cast<FT>(CGAL::sqrt(CGAL::to_double(squared_distance(corner, points[segment_index][i]))));
				}

				int closest_index = -1; FT min_dist = FT(1000000000);
				for (size_t i = 0; i < dist.size(); ++i) {

					if (dist[i] < min_dist) {

						closest_index = static_cast<int>(i);
						min_dist = dist[i];
					}
				}

				assert(closest_index != -1);
				return closest_index;
			}

			template<class ContainerT, class Key>
			bool is_valid_key(const ContainerT &ct, const Key &key) const {
				return ct.count(key);
			}

			template<class Points>
			int compute_delaunay(const Points &points, CDT &cdt) const {

				Log log;

				auto number_of_faces = -1;
				assert(!points.empty());
				
				cdt.clear();

				using Point_iterator = typename Points::const_iterator;


				// Insert all points.
				for (Point_iterator pit = points.begin(); pit != points.end(); ++pit) {
					
					const auto &point = (*pit).second;
					cdt.insert(point);
				}


				// Create CDT.
				// CGAL::make_conforming_Delaunay_2(cdt);
				number_of_faces = cdt.number_of_faces();


				// Save CDT.
				log.save_cdt_obj(cdt, "tmp" + std::string(PSR) + "cdt");

				// log.save_cdt_ply(cdt, "tmp" + std::string(PSR) + "visibility_before", "in");
				// log.save_cdt_ply(cdt, "tmp" + std::string(PSR) + "buildings_before" , "bu");

				return number_of_faces;
			}

			void insert_constraint(CDT &cdt, const Edge_iterator &edge_handle) const {

				const Face_handle face = edge_handle->first;
				const int vertex_index = edge_handle->second;

				const Vertex_handle &source = face->vertex(cdt.ccw(vertex_index));
				const Vertex_handle &target = face->vertex( cdt.cw(vertex_index));

				cdt.insert_constraint(source, target);
			}

			bool is_boundary_edge(const CDT &cdt, const Edge_iterator &edge_handle) const {

				const int vertex_index = edge_handle->second;

				const Face_handle face_1 = edge_handle->first;
				const Face_handle face_2 = face_1->neighbor(vertex_index);

				if (cdt.is_infinite(face_1) || cdt.is_infinite(face_2)) return true;
				return false;
			}


			//////////////////////////
			// Bounding box functions!

			template<class Indices, class Box>
			void return_bounding_box(const Plane_3 &plane, const Input &input, const Indices &indices, Box &box) {

				const FT big_value = FT(100000000000000);

				FT minx = big_value, miny = big_value, maxx = -big_value, maxy = -big_value, z;
				for (size_t i = 0; i < indices.size(); ++i) {
					
					const Point_3 &p = input.point(indices[i]);
					const Point_3 projected = plane.projection(p);

					minx = CGAL::min(minx, projected.x());
					miny = CGAL::min(miny, projected.y());

					maxx = CGAL::max(maxx, projected.x());
					maxy = CGAL::max(maxy, projected.y());

					z = projected.z();
				}

				box.bbmin = Point_3(minx, miny, z);
				box.bbmax = Point_3(maxx, maxy, z);
			}

			void add_bbox_to(const CDT &cdt, const Input &input, CDT &cdt_with_bbox) {

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

			template<class Ground, class Ground_point>
			void compute_ground_bbox(const Input &input, Ground &ground_bbox) const {

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

			void compute_bounding_box(const Input &input, std::vector<Point_2> &bbox) const {

				bbox.clear();
				bbox.resize(4);

				const FT big_value = FT(1000000000); // change it

				FT minx =  big_value, miny =  big_value;
				FT maxx = -big_value, maxy = -big_value;

				for (typename Input::const_iterator it = input.begin(); it != input.end(); ++it) {
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


			//////////////////////////////////
			// Ground plane fitting functions!

			// Not efficient since I need to copy all ground points.
			template<class Indices>
			void fit_ground_plane(const Input &input, const Indices &ground_idxs, Plane_3 &ground_plane) const {

				using Local_Kernel = CGAL::Simple_cartesian<double>;
				using Point_3ft    = Local_Kernel::Point_3;
				using Plane_3ft    = Local_Kernel::Plane_3;

				std::vector<Point_3ft> tmp_ground(ground_idxs.size());
				for (size_t i = 0; i < ground_idxs.size(); ++i) {
				
					const Point_3 &p = input.point(ground_idxs[i]);

					const double x = CGAL::to_double(p.x());
					const double y = CGAL::to_double(p.y());
					const double z = CGAL::to_double(p.z());

					tmp_ground[i] = Point_3ft(x, y, z);
				}

				Plane_3ft tmp_plane;
				CGAL::linear_least_squares_fitting_3(tmp_ground.begin(), tmp_ground.end(), tmp_plane, CGAL::Dimension_tag<0>());
				ground_plane = Plane_3(static_cast<FT>(tmp_plane.a()), static_cast<FT>(tmp_plane.b()), static_cast<FT>(tmp_plane.c()), static_cast<FT>(tmp_plane.d()));
			}


			/////////////////////
			// Mapping functions!

			template<class Container_2D, class Face_points_map>
			int get_2d_input_and_face_points_map(const CDT &cdt, const Input &input_3d, Container_2D &input_2d, Face_points_map &fp_map, const bool silent) const {

				Log log_all, log_in_cdt;

				typename CDT::Locate_type locate_type;
				int locate_index_stub = -1;

				Label_map class_labels;
				boost::tie(class_labels, boost::tuples::ignore) = input_3d.template property_map<Label>("label");

				input_2d.clear();
				input_2d.resize(input_3d.number_of_points());

				size_t point_index = 0;
				for (typename Input::const_iterator it = input_3d.begin(); it != input_3d.end(); ++it) {

					const Point_index pi = *it;

					const Point_3 &p = input_3d.point(pi);
					const Point_2 &q = Point_2(p.x(), p.y());

					input_2d[point_index] = std::make_pair(q, class_labels[pi]);
					log_all.out << q << " " << 0 << std::endl;

					const Face_handle face_handle = cdt.locate(q, locate_type, locate_index_stub);
					if (locate_type == CDT::FACE || locate_type == CDT::EDGE || locate_type == CDT::VERTEX) {
						
						fp_map[face_handle].push_back(pi);

						log_in_cdt.out << q << " " << 0 << std::endl;
						++point_index;
					}
				}

				if (!silent) log_all.save("tmp" + std::string(PSR) + "lod_0_1" + std::string(PSR) + "input_2d_all", ".xyz");
				
				// log_in_cdt.save("tmp" + std::string(PSR) + "input_2d_in_cdt", ".xyz");
				
				return point_index;
			}


			//////////////////////////////////
			// Segment-line related functions!

			// Not efficient since I need to copy all ground points.
			template<class Projected_points, class Planes, class Lines>
			int fit_lines_to_projected_points(const Projected_points &points, const Planes &planes, Lines &lines) const {

				using Plane_iterator = typename Planes::const_iterator;

      			using Local_Kernel = CGAL::Simple_cartesian<double>;
				using Point_2ft    = Local_Kernel::Point_2;
				using Line_2ft     = Local_Kernel::Line_2;

				auto number_of_fitted_lines = 0;
				std::vector<Point_2ft> tmp_points;

				lines.clear();
				lines.resize(planes.size());

				for (Plane_iterator it = planes.begin(); it != planes.end(); ++it, ++number_of_fitted_lines) {
					const auto num_points = (*it).second.size();

					tmp_points.clear();
					tmp_points.resize(num_points);

					for (size_t i = 0; i < num_points; ++i) {
						
						const auto index = (*it).second[i];
						const Point_2 &p = points.at(index);

						const double x = CGAL::to_double(p.x());
						const double y = CGAL::to_double(p.y());

						tmp_points[i] = Point_2ft(x, y);
					}

					Line_2ft tmp_line;
					CGAL::linear_least_squares_fitting_2(tmp_points.begin(), tmp_points.end(), tmp_line, CGAL::Dimension_tag<0>());
					lines[number_of_fitted_lines] = Line_2(static_cast<FT>(tmp_line.a()), static_cast<FT>(tmp_line.b()), static_cast<FT>(tmp_line.c()));
				}
				return number_of_fitted_lines;
			}

			// It may have precision problems: e.g. Line = (1 1.0e+24) -- (1 1.0e+24) or nan!
			template<class Projected_points, class Planes, class Lines, class Segments>
			int create_segments_from_lines(const Projected_points &points, const Planes &planes, const Lines &lines, Segments &segments) const {

				using Plane_iterator = typename Planes::const_iterator;

				segments.clear();
				segments.resize(lines.size());

				std::cout.precision(20);

				auto number_of_segments = 0;
				for (Plane_iterator it = planes.begin(); it != planes.end(); ++it, ++number_of_segments) {

					const auto num_points = (*it).second.size();
					const auto big = +FT(100000000000000);

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
					const auto eps = FT(1) / FT(1000000);
					if ((v1.y() < FT(0) && v2.y() >= FT(0) && CGAL::abs(v1.y() - v2.y()) > eps) ||
						(v2.y() < FT(0) && v1.y() >= FT(0) && CGAL::abs(v1.y() - v2.y()) > eps) ||
						(v1.x() < FT(0) && v2.x() >= FT(0) && CGAL::abs(v1.x() - v2.x()) > eps) ||
						(v2.x() < FT(0) && v1.x() >= FT(0) && CGAL::abs(v1.x() - v2.x()) > eps)) {

						segments[number_of_segments] = Segment_2(Point_2(minx, maxy), Point_2(maxx, miny));
					}
				}

				assert(static_cast<size_t>(number_of_segments) == lines.size());
				return number_of_segments;
			}

			// My custom function to handle precision problems when projecting points.
			inline Point_2 project(const Line_2 &line, const Point_2 &p) const {
				return m_simple_utils.project_onto_line(line, p);
			}

			inline Point_2 project(const Segment_2 &segment, const Point_2 &p) const {
				return m_simple_utils.project_onto_line(segment, p);
			}
		};
	}
}

#endif // CGAL_LEVEL_OF_DETAIL_UTILS_H