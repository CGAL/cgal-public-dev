#ifndef CGAL_LEVEL_OF_DETAIL_DISTORTION_H
#define CGAL_LEVEL_OF_DETAIL_DISTORTION_H

#if defined(WIN32) || defined(_WIN32) 
#define PS "\\"
#else 
#define PS "/"
#endif

// STL includes.
#include <map>
#include <string>
#include <vector>
#include <cassert>
#include <iostream>

#define BOOST_VERSION 104500

// CGAL includes.
#include <CGAL/number_utils.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_tree.h>

// New CGAL includes.
#include <CGAL/Mylog/Mylog.h>
#include <CGAL/Selector/Level_of_detail_selector.h>
#include <CGAL/Selector/Level_of_detail_selection_strategy.h>
#include <CGAL/Utils/Level_of_detail_utils_simple.h>
#include <CGAL/Level_of_detail_enum.h>

namespace CGAL {

	namespace LOD {

		template<class KernelTraits, class InputContainer, class LodReconstruction>
		class Level_of_detail_distortion {

		public:
			typedef KernelTraits   	  Kernel;
			typedef InputContainer 	  Container;
			typedef LodReconstruction LODS;

			typedef typename LodReconstruction::Mesh 			  Mesh;
			typedef typename LodReconstruction::Mesh_facet_colors Mesh_facet_colors;

			using Vertex_handle   = typename Mesh::Vertex_handle;
			using Halfedge_handle = typename Mesh::Halfedge_const_handle;
			using Facet_handle    = typename Mesh::Facet_const_handle;
			using Faces 	      = std::vector< std::vector<Facet_handle> >;
			
			typedef typename Kernel::FT 	    FT;
			typedef typename Kernel::Point_3    Point_3;
			typedef typename Kernel::Vector_3   Vector_3;
			typedef typename Kernel::Plane_3    Plane_3;
			typedef typename Kernel::Point_2    Point_2;
			typedef typename Kernel::Triangle_2 Triangle_2;

			using Index   = int;
			using Indices = std::vector<Index>;

			typedef CGAL::LOD::Level_of_detail_building_interior<Kernel, Container> Roofs_points_selection_strategy;
			typedef CGAL::LOD::Level_of_detail_building_boundary<Kernel, Container> Walls_points_selection_strategy;
			
			typedef CGAL::LOD::Level_of_detail_selector<Kernel, Roofs_points_selection_strategy> Roofs_points_selector;
			typedef CGAL::LOD::Level_of_detail_selector<Kernel, Walls_points_selection_strategy> Walls_points_selector;

			typedef CGAL::AABB_face_graph_triangle_primitive<Mesh> AB_primitive;
			typedef CGAL::AABB_traits<Kernel, AB_primitive> 	   AB_traits;
			typedef CGAL::AABB_tree<AB_traits> 					   AB_tree;

			typedef typename AB_tree::Point_and_primitive_id AB_point_and_primitive_id;
			typename Kernel::Compute_squared_distance_3 squared_distance;
			using Log = CGAL::LOD::Mylog;

			typedef CGAL::LOD::Level_of_detail_utils_simple<Kernel> Simple_utils;
			struct Bounding_box { Point_3 bbmin, bbmax; };

			using Bounding_boxes = std::vector<Bounding_box>;

			Level_of_detail_distortion(const Container &input, const LODS &lods) : 
			m_input(input), m_lods(lods),
			m_distortion(-FT(1)), m_num_roofs_points(-FT(1)), m_num_walls_points(-FT(1)),
			m_debug(false) { }

			void estimate() {

				Mesh roofs_mesh, walls_mesh;
				get_roofs(roofs_mesh);
				get_walls(walls_mesh);

				Indices roofs_point_indices, walls_point_indices;
				get_roofs_points(roofs_point_indices);
				get_walls_points(walls_point_indices);

				std::vector<Point_3> roofs_translated_points;
				get_translated_points(roofs_mesh, roofs_point_indices, roofs_translated_points, "roofs_translated_points_for_distortion");
				m_num_roofs_points = roofs_translated_points.size();

				std::vector<Point_3> walls_translated_points;
				get_translated_points(walls_mesh, walls_point_indices, walls_translated_points, "walls_translated_points_for_distortion");
				m_num_walls_points = walls_translated_points.size();

				compute_metrics(roofs_mesh, roofs_translated_points, walls_mesh, walls_translated_points);
			}

			FT get(const Quality_data_type type = Quality_data_type::DST_AVG) {

				m_distortion = fit_data(type);

				assert(m_distortion >= FT(0));
				return m_distortion;
			}

			const std::vector<FT> &get_roofs_metrics() const {

				assert(!m_roofs_metrics.empty());
				return m_roofs_metrics;
			}

			const std::vector<FT> &get_walls_metrics() const {
				
				assert(!m_walls_metrics.empty());
				return m_walls_metrics;
			}

			FT number_of_roofs_points() const {
				
				assert(m_num_roofs_points >= FT(0));
				return m_num_roofs_points;
			}

			FT number_of_walls_points() const {
				
				assert(m_num_walls_points >= FT(0));
				return m_num_walls_points;
			}

		private:
			const Container &m_input;
			const LODS 		&m_lods;

			FT   m_distortion, m_num_roofs_points, m_num_walls_points;
			bool m_debug;

			std::vector<FT> m_roofs_metrics;
			std::vector<FT> m_walls_metrics;
			
			Simple_utils m_simple_utils;

			void get_roofs(Mesh &roofs) {
				m_lods.get_roofs(roofs);

				if (m_debug) {
					Log roofs_saver; Mesh_facet_colors stub;
					roofs_saver.save_mesh_as_ply(roofs, stub, "tmp" + std::string(PS) + "roofs_for_distortion", false);
				}
			}

			void get_walls(Mesh &walls) {
				m_lods.get_walls(walls);

				if (m_debug) {
					Log walls_saver; Mesh_facet_colors stub;
					walls_saver.save_mesh_as_ply(walls, stub, "tmp" + std::string(PS) + "walls_for_distortion", false);
				}
			}

			void get_roofs_points(Indices &roofs_point_indices) {
				roofs_point_indices.clear();
				
				Roofs_points_selector selector;
				selector.select_elements(m_input, std::back_inserter(roofs_point_indices));
	
				if (m_debug) {
					Log roofs_points_saver;
					roofs_points_saver.export_points_using_indices(m_input, roofs_point_indices, "tmp" + std::string(PS) + "roofs_points_for_distortion");
				}
			}

			void get_walls_points(Indices &walls_point_indices) {
				walls_point_indices.clear();
				
				Walls_points_selector selector;
				selector.select_elements(m_input, std::back_inserter(walls_point_indices));
	
				if (m_debug) {
					Log walls_points_saver;
					walls_points_saver.export_points_using_indices(m_input, walls_point_indices, "tmp" + std::string(PS) + "walls_points_for_distortion");
				}
			}

			void get_translated_points(const Mesh &mesh, const Indices &point_indices, std::vector<Point_3> &translated_points, const std::string &name) {

				const Vector_3 translation = get_translation(mesh, point_indices);
				translate_points(translation, point_indices, translated_points, name);
			}

			Vector_3 get_translation(const Mesh &mesh, const Indices &point_indices) {

				const Point_3 mesh_barycentre   = get_mesh_barycentre(mesh);
				const Point_3 points_barycentre = get_points_barycentre(point_indices);

				return Vector_3(mesh_barycentre, points_barycentre);
			}

			Point_3 get_mesh_barycentre(const Mesh &mesh) {

				FT num_points = FT(0);
				FT x = FT(0), y = FT(0), z = FT(0);

				for (typename Mesh::Vertex_const_iterator vit = mesh.vertices_begin(); vit != mesh.vertices_end(); ++vit, num_points += FT(1)) {
					const Point_3 &p = vit->point();

					x += p.x();
					y += p.y();
					z += p.z();
				}
				x /= num_points;
				y /= num_points;
				z /= num_points;

				return Point_3(x, y, z);
			}

			Point_3 get_points_barycentre(const Indices &point_indices) {

				FT num_points = FT(0);
				FT x = FT(0), y = FT(0), z = FT(0);

				for (size_t i = 0; i < point_indices.size(); ++i, num_points += FT(1)) {
					const Point_3 &p = m_input.point(point_indices[i]);

					x += p.x();
					y += p.y();
					z += p.z();
				}
				x /= num_points;
				y /= num_points;
				z /= num_points;

				return Point_3(x, y, z);
			}

			void translate_points(const Vector_3 &translation, const Indices &point_indices, std::vector<Point_3> &translated_points, const std::string &name) {

				translated_points.clear();
				translated_points.resize(point_indices.size());

				for (size_t i = 0; i < point_indices.size(); ++i) {
					const Point_3 &p = m_input.point(point_indices[i]);
					const FT z = p.z() - translation.z();

					translated_points[i] = Point_3(p.x(), p.y(), z);
				}

				if (m_debug) {
					Log saver;
					saver.export_points(translated_points, "tmp" + std::string(PS) + name);
				}
			}

			void compute_metrics(
				const Mesh &roofs_mesh, const std::vector<Point_3> &roofs_translated_points,
				const Mesh &walls_mesh, const std::vector<Point_3> &walls_translated_points) {

				compute_l2_distances(roofs_mesh, roofs_translated_points, walls_mesh, walls_translated_points);
			}

			void compute_l2_distances(
				const Mesh &roofs_mesh, const std::vector<Point_3> &roofs_translated_points,
				const Mesh &walls_mesh, const std::vector<Point_3> &walls_translated_points) {

				compute_mesh_points_l2_distances_to_roofs(roofs_mesh, roofs_translated_points);
				compute_mesh_points_l2_distances_to_walls(walls_mesh, walls_translated_points);
			}

			void compute_mesh_points_l2_distances_to_roofs(const Mesh &mesh, const std::vector<Point_3> &points) {

				Faces roofs_faces;
				m_lods.get_roofs_faces(roofs_faces);

				Bounding_boxes boxes;
				compute_bounding_boxes(roofs_faces, boxes);

				std::vector<Point_3> updated_points;
				for (size_t i = 0; i < points.size(); ++i) {
					const Point_3 &query = points[i];
					
					std::vector<int> box_indices;
					find_bounding_boxes(query, boxes, box_indices);
					
					if (box_indices.empty()) {
						updated_points.push_back(query); continue;
					}

					bool found = false;
					for (size_t j = 0; j < box_indices.size(); ++j) {

						const Plane_3 plane = get_plane_from_bounding_box(boxes[box_indices[j]].bbmin, boxes[box_indices[j]].bbmax);
						const Point_3 projected = plane.projection(query);

						const int face_index = find_face(projected, roofs_faces[box_indices[j]]);
						if (face_index >= 0) {
							found = true; break;
						}
					}

					if (!found) {
						updated_points.push_back(query); continue;
					}
				}

				if (m_debug) {
					std::cout << std::endl << "number of original - updated points: " << points.size() << " - " << updated_points.size() << std::endl << std::endl;
					
					Log saver;
					saver.export_points(updated_points, "tmp" + std::string(PS) + "outlier_points_for_roofs_in_distortion");
				}
				compute_updated_mesh_points_l2_distances_to_roofs(mesh, updated_points);
			}

			void compute_bounding_boxes(const Faces &faces, Bounding_boxes &boxes) {

				assert(!faces.empty());
				boxes.resize(faces.size());

				std::vector<Point_3> points;
				for (size_t i = 0; i < faces.size(); ++i) {
					
					get_points(faces[i], points);
					m_simple_utils.compute_bounding_box_in_3d(boxes[i].bbmin, boxes[i].bbmax, points);

					assert(boxes[i].bbmin.z() == boxes[i].bbmax.z());
				}

				if (m_debug) {
					Log log;
					log.save_bounding_boxes_as_ply<Bounding_boxes, Point_3>(boxes, "tmp" + std::string(PS) + "boxes");
				}
			}

			void get_points(const std::vector<Facet_handle> &faces, std::vector<Point_3> &points) {
				assert(!faces.empty());

				points.clear();
				for (size_t i = 0; i < faces.size(); ++i) {
					
					Halfedge_handle he = faces[i]->halfedge();
					points.push_back(he->vertex()->point());

					he = he->next();
					points.push_back(he->vertex()->point());

					he = he->next();
					points.push_back(he->vertex()->point());
				}
			}

			void find_bounding_boxes(const Point_3 &query, const Bounding_boxes &boxes, std::vector<int> &box_indices) {
				
				assert(!boxes.empty());
				assert(box_indices.empty());

				for (size_t i = 0; i < boxes.size(); ++i)
					if (is_inside_box(query, boxes[i].bbmin, boxes[i].bbmax)) 
						box_indices.push_back(static_cast<int>(i));
			}

			bool is_inside_box(const Point_3 &query, const Point_3 &minp, const Point_3 &maxp) {

				if (query.x() > minp.x() && query.x() < maxp.x() && 
					query.y() > minp.y() && query.y() < maxp.y()) return true;

				return false;
			}

			Plane_3 get_plane_from_bounding_box(const Point_3 &minp, const Point_3 &maxp) {
				assert(minp.z() == maxp.z());

				const Point_3 a = minp;
				const Point_3 b = Point_3(maxp.x(), minp.y(), minp.z());
				const Point_3 c = maxp;

				return Plane_3(a, b, c);
			}

			int find_face(const Point_3 &query, const std::vector<Facet_handle> &faces) {
				assert(!faces.empty());

				for (size_t i = 0; i < faces.size(); ++i)
					if (belongs_to_face(query, faces[i]))
						return static_cast<int>(i);

				return -1;
			}

			bool belongs_to_face(const Point_3 &query, const Facet_handle &fh) {

				std::vector<Point_3> points;
				get_points_from_face_handle(fh, points);

				assert(points.size() == 3);
				assert(points[0].z() == points[1].z() && points[1].z() == points[2].z());

				const Triangle_2 triangle = Triangle_2(
					Point_2(points[0].x(), points[0].y()),
					Point_2(points[1].x(), points[1].y()),
					Point_2(points[2].x(), points[2].y())
					);

				const Point_2 new_query = Point_2(query.x(), query.y());

				if (triangle.has_on_bounded_side(new_query) || triangle.has_on_boundary(new_query)) return true;
				return false;
			}

			void get_points_from_face_handle(const Facet_handle &fh, std::vector<Point_3> &points) {
				
				points.clear();
				Halfedge_handle he = fh->halfedge();

				const Point_3 p1 = he->vertex()->point();

				he = he->next();
				const Point_3 p2 = he->vertex()->point();

				he = he->next();
				const Point_3 p3 = he->vertex()->point();

				points.push_back(p1);
				points.push_back(p2);
				points.push_back(p3);
			}

			void compute_updated_mesh_points_l2_distances_to_roofs(const Mesh &mesh, const std::vector<Point_3> &points) {
				
				assert(m_num_roofs_points > FT(0));
				assert(points.size() <= m_num_roofs_points);

				m_roofs_metrics.clear();
				m_roofs_metrics.resize(static_cast<size_t>(CGAL::to_double(m_num_roofs_points)));
				for (size_t i = 0; i < m_roofs_metrics.size(); ++i) m_roofs_metrics[i] = FT(0);

				Mesh tmp_mesh = mesh;
				for (Vertex_handle vh = tmp_mesh.vertices_begin(); vh != tmp_mesh.vertices_end(); ++vh) {

					Point_3 &p = vh->point();
					p = Point_3(p.x(), p.y(), FT(0));
				}

				AB_tree aabb_tree(faces(tmp_mesh).first, faces(tmp_mesh).second, tmp_mesh);

				Point_3 closest_point;
				for (size_t i = 0; i < points.size(); ++i) {
					
					const Point_3 query = Point_3(points[i].x(), points[i].y(), FT(0));
					closest_point = aabb_tree.closest_point(query);

					const FT squared_dist = squared_distance(query, closest_point);
					const FT distance = static_cast<FT>(CGAL::sqrt(CGAL::to_double(squared_dist)));

					m_roofs_metrics[i] = distance;
				}
			}

			void compute_mesh_points_l2_distances_to_walls(const Mesh &mesh, const std::vector<Point_3> &points) {

				assert(m_num_walls_points > FT(0));
				assert(points.size() <= m_num_walls_points);

				m_walls_metrics.clear();
				m_walls_metrics.resize(static_cast<size_t>(CGAL::to_double(m_num_walls_points)));
				
				for (size_t i = 0; i < m_walls_metrics.size(); ++i) m_walls_metrics[i] = FT(0);
				assert(m_walls_metrics.size() == points.size());

				Mesh tmp_mesh = mesh;
				for (Vertex_handle vh = tmp_mesh.vertices_begin(); vh != tmp_mesh.vertices_end(); ++vh) {

					Point_3 &p = vh->point();
					p = Point_3(p.x(), p.y(), FT(0));
				}

				AB_tree aabb_tree(faces(tmp_mesh).first, faces(tmp_mesh).second, tmp_mesh);

				Point_3 closest_point;
				for (size_t i = 0; i < points.size(); ++i) {
					
					const Point_3 query = Point_3(points[i].x(), points[i].y(), FT(0));
					closest_point = aabb_tree.closest_point(query);

					const FT squared_dist = squared_distance(query, closest_point);
					const FT distance = static_cast<FT>(CGAL::sqrt(CGAL::to_double(squared_dist)));

					m_walls_metrics[i] = distance;
				}
			}

			Plane_3 get_plane_from_face_handle(const Facet_handle &fh) {
				
				std::vector<Point_3> points;
				get_points_from_face_handle(fh, points);

				assert(points.size() == 3);
				return Plane_3(points[0], points[1], points[2]);
			}

			FT fit_data(const Quality_data_type type) const {

				switch (type) {
					case Quality_data_type::DST_MIN:
						return get_min_from_data();

					case Quality_data_type::DST_AVG:
						return get_avg_from_data();

					case Quality_data_type::DST_MAX:
						return get_max_from_data();

					default:
						assert(!"Wrong fitting type!");
						return -FT(1);
				}
			}

			FT get_min_from_data() const {

				const FT roofs_mind = get_min(m_roofs_metrics);
				const FT walls_mind = get_min(m_walls_metrics);

				const FT mind = CGAL::min(roofs_mind, walls_mind);
				return mind;
			}

			FT get_avg_from_data() const {

				assert(m_num_roofs_points > FT(0));
				assert(m_num_walls_points > FT(0));

				const FT roofs_avg = get_sum(m_roofs_metrics) / FT(m_num_roofs_points);
				const FT walls_avg = get_sum(m_walls_metrics) / FT(m_num_walls_points);

				if (m_debug) std::cout << "roofs avg " << roofs_avg << ", walls avg " << walls_avg << std::endl << std::endl;

				const FT avg = (roofs_avg + walls_avg) / FT(2);
				return avg;
			}

			FT get_max_from_data() const {

				const FT roofs_maxd = get_max(m_roofs_metrics);
				const FT walls_maxd = get_max(m_walls_metrics);

				const FT maxd = CGAL::max(roofs_maxd, walls_maxd);
				return maxd;
			}

			FT get_min(const std::vector<FT> &data) const {
				assert(!data.empty());

				FT mind = FT(1000000000000);
				for (size_t i = 0; i < data.size(); ++i) mind = CGAL::min(mind, data[i]);

				return mind;
			}

			FT get_sum(const std::vector<FT> &data) const {
				assert(!data.empty());

				FT sum = FT(0);
				for (size_t i = 0; i < data.size(); ++i) sum += data[i];

				return sum;
			}

			FT get_max(const std::vector<FT> &data) const {
				assert(!data.empty());

				FT maxd = -FT(1000000000000);
				for (size_t i = 0; i < data.size(); ++i) maxd = CGAL::max(maxd, data[i]);

				return maxd;
			}
		};
	}
}

#endif // CGAL_LEVEL_OF_DETAIL_DISTORTION_H