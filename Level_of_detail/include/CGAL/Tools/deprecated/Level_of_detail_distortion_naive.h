#ifndef CGAL_LEVEL_OF_DETAIL_DISTORTION_NAIVE_H
#define CGAL_LEVEL_OF_DETAIL_DISTORTION_NAIVE_H

#if defined(WIN32) || defined(_WIN32) 
#define PS "\\"
#else 
#define PS "/"
#endif

// STL includes
#include <string>
#include <vector>
#include <cassert>
#include <iostream>

// CGAL includes.
#include <CGAL/number_utils.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_tree.h>

// New CGAL includes.
#include <CGAL/Mylog/Mylog.h>
#include <CGAL/Selector/Level_of_detail_selector.h>
#include <CGAL/Selector/Level_of_detail_selection_strategy.h>

namespace CGAL {

	namespace LOD {

		template<class KernelTraits, class InputContainer, class LodReconstruction>
		class Level_of_detail_distortion_naive {

		public:
			typedef KernelTraits   	  Kernel;
			typedef InputContainer 	  Container;
			typedef LodReconstruction LODS;

			typedef typename LodReconstruction::Mesh Mesh;
			typedef typename LodReconstruction::Mesh_facet_colors Mesh_facet_colors;
			
			typedef typename Kernel::FT 	  FT;
			typedef typename Kernel::Point_3  Point_3;
			typedef typename Kernel::Vector_3 Vector_3;

			using Index   = int;
			using Indices = std::vector<Index>;

			typedef CGAL::LOD::Level_of_detail_building_boundary<Kernel, Container> Boundary_points_strategy;
			typedef CGAL::LOD::Level_of_detail_building_interior<Kernel, Container> Interior_points_strategy;

			typedef CGAL::LOD::Level_of_detail_selector<Kernel, Boundary_points_strategy> Boundary_points_selector;
			typedef CGAL::LOD::Level_of_detail_selector<Kernel, Interior_points_strategy> Interior_points_selector;

			typedef CGAL::AABB_face_graph_triangle_primitive<Mesh> AB_primitive;
			typedef CGAL::AABB_traits<Kernel, AB_primitive> 	   AB_traits;
			typedef CGAL::AABB_tree<AB_traits> 					   AB_tree;

			typename Kernel::Compute_squared_distance_3 squared_distance;
			using Log = CGAL::LOD::Mylog;


			Level_of_detail_distortion_naive(const Container &input, const LODS &lods) : 
			m_input(input), m_lods(lods), 
			m_distortion(-FT(1)), m_save_info(false) { }


			void estimate() {
				
				Mesh roofs_mesh, walls_mesh;
				get_roofs(roofs_mesh);
				get_walls(walls_mesh);

				Indices roofs_point_indices, walls_point_indices;
				get_roofs_points(roofs_point_indices);
				get_walls_points(walls_point_indices);

				std::vector<Point_3> roofs_translated_points;
				get_translated_points(roofs_mesh, roofs_point_indices, roofs_translated_points, "roofs_translated_points");

				std::vector<Point_3> walls_translated_points;
				get_translated_points(walls_mesh, walls_point_indices, walls_translated_points, "walls_translated_points");

				const FT average_distance_to_roofs = compute_mesh_points_average_distance(roofs_mesh, roofs_translated_points);
				const FT average_distance_to_walls = compute_mesh_points_average_distance(walls_mesh, walls_translated_points);

				m_distortion = (average_distance_to_roofs + average_distance_to_walls) / FT(2);
			}

			FT get() const {
				return m_distortion;
			}

		private:
			const Container &m_input;
			const LODS 		&m_lods;

			FT m_distortion;
			bool m_save_info;


			void get_roofs(Mesh &roofs) {
				m_lods.get_roofs(roofs);

				if (m_save_info) {
					Log roofs_saver; Mesh_facet_colors stub;
					roofs_saver.save_mesh_as_ply(roofs, stub, "tmp" + std::string(PS) + "roofs_for_complexity", false);
				}
			}

			void get_walls(Mesh &walls) {
				m_lods.get_walls(walls);

				if (m_save_info) {
					Log walls_saver; Mesh_facet_colors stub;
					walls_saver.save_mesh_as_ply(walls, stub, "tmp" + std::string(PS) + "walls_for_complexity", false);
				}
			}

			void get_roofs_points(Indices &roofs_points) {
				roofs_points.clear();
				
				Interior_points_selector selector;
				selector.select_elements(m_input, std::back_inserter(roofs_points));
	
				if (m_save_info) {
					Log roofs_points_saver;
					roofs_points_saver.export_points_using_indices(m_input, roofs_points, "tmp" + std::string(PS) + "roofs_points_for_complexity");
				}
			}

			void get_walls_points(Indices &walls_points) {
				walls_points.clear();
				
				Boundary_points_selector selector;
				selector.select_elements(m_input, std::back_inserter(walls_points));
	
				if (m_save_info) {
					Log walls_points_saver;
					walls_points_saver.export_points_using_indices(m_input, walls_points, "tmp" + std::string(PS) + "walls_points_for_complexity");
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

				if (m_save_info) {
					Log saver;
					saver.export_points(translated_points, "tmp" + std::string(PS) + name);
				}
			}

			FT compute_mesh_points_average_distance(const Mesh &mesh, const std::vector<Point_3> &points) {

				AB_tree aabb_tree(faces(mesh).first, faces(mesh).second, mesh);
				Point_3 closest_point;

				FT average_distance = FT(0);
				FT num_points = FT(0);

				for (size_t i = 0; i < points.size(); ++i, num_points += FT(1)) {
					
					const Point_3 &query = points[i];
					closest_point = aabb_tree.closest_point(query);

					const FT squared_dist = squared_distance(query, closest_point);
					const FT distance = static_cast<FT>(CGAL::sqrt(CGAL::to_double(squared_dist)));

					average_distance += distance;
				}

				average_distance /= num_points;
				return average_distance;
			}
		};
	}
}

#endif // CGAL_LEVEL_OF_DETAIL_DISTORTION_NAIVE_H