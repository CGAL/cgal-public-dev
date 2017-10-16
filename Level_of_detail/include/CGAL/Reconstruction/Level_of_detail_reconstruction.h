#ifndef CGAL_LEVEL_OF_DETAIL_RECONSTRUCTION_H
#define CGAL_LEVEL_OF_DETAIL_RECONSTRUCTION_H

// STL includes.
#include <vector>
#include <iostream>
#include <cassert>
#include <string>
#include <map>

// CGAL includes.
#include <CGAL/number_utils.h>
#include <CGAL/utils_classes.h>
#include <CGAL/utils.h>
#include <CGAL/IO/Color.h>
#include <CGAL/Unique_hash_map.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/Random.h>

// New CGAL includes.
#include <CGAL/Mylog/Mylog.h>
#include <CGAL/Level_of_detail_enum.h>

namespace CGAL {

	namespace LOD {

		// Mesh builder for the LOD1 reconstruction below.
		template<class KernelTraits, class HDSInput, class CDTInput, class BuildingsInput, class ColorFacetHandle>
		class Build_mesh : public Modifier_base<HDSInput> {
		
		public:
			typedef KernelTraits 	 Kernel;
    		typedef HDSInput 		 HDS;
    		typedef CDTInput 		 CDT;
    		typedef BuildingsInput 	 Buildings;
    		typedef ColorFacetHandle Color_facet_handle;

    		typedef typename HDS::Vertex   Vertex;
    		typedef typename Vertex::Point Point;

    		enum class Build_type { TEST, CDT_AND_BUILDINGS };

    		using FT = typename Kernel::FT;

			using Color 	   = CGAL::Color;
			using Facet_colors = std::map<Color_facet_handle, Color>;
			using Builder 	   = CGAL::Polyhedron_incremental_builder_3<HDS>;
    		
    		using Log = CGAL::LOD::Mylog;
    		using Building_iterator = typename Buildings::const_iterator;
    		using Ground = std::vector<Point>;

    		enum class Builder_type { LOD0, LOD1 };

    		Build_mesh(const CDT &cdt, const Buildings &buildings, Facet_colors &facet_colors) : 
    		m_build_type(Build_type::CDT_AND_BUILDINGS), 
    		m_cdt(cdt), 
    		m_buildings(buildings), 
    		m_facet_colors(facet_colors), 
    		m_index_counter(0), 
    		m_ground_set(false), 
    		m_builder_type(Builder_type::LOD1),
    		m_use_boundaries(true) { }

			void operator()(HDS &hds) {

				Builder builder(hds, false);

				switch (m_builder_type) {

					case Builder_type::LOD0:
						build_lod0(builder);
						break;

					case Builder_type::LOD1:
						build_lod1(builder);
						break;

					default:
						assert(!"Wrong builder type!");
						break;
				}
			}

			void set_builder_type(const Builder_type builder_type) {
				m_builder_type = builder_type;
			}

			void set_ground(const Ground &ground) {
				m_ground = ground;
				m_ground_set = true;
			}

			void use_boundaries(const bool new_state) {
				m_use_boundaries = new_state;
			}

		private:
			const Build_type m_build_type;
			const CDT 		 &m_cdt;
			const Buildings  &m_buildings;
			
			Facet_colors &m_facet_colors;
			CGAL::Random m_rand;
			
			size_t m_index_counter;
			Ground m_ground;
			
			bool m_ground_set;
			Builder_type m_builder_type;
			bool m_use_boundaries;
			
			void build_lod0(Builder &builder) {
				assert(m_build_type == Build_type::CDT_AND_BUILDINGS);

				const size_t expected_num_vertices = estimate_number_of_vertices_lod0();
				const size_t expected_num_faces    = estimate_number_of_faces_lod0();


				// Start surface building.
				m_index_counter = 0;
				builder.begin_surface(expected_num_vertices, expected_num_faces);


				// Add all buildings.
				for (Building_iterator bit = m_buildings.begin(); bit != m_buildings.end(); ++bit) {
				
					const auto &building = (*bit).second;
					add_new_building_lod0(building, builder);
				}


				// Add ground.
				assert(m_ground_set);
				add_ground(builder);


				// End surface building.
				builder.end_surface();
			}

			void build_lod1(Builder &builder) {

				switch (m_build_type) {

					case Build_type::TEST:
						build_test_data_lod1(builder);
						break;

					case Build_type::CDT_AND_BUILDINGS:
						build_from_cdt_and_buildings_lod1(builder);
						break;

					default:
						assert(!"Wrong build type!");
						break;
				}
			}

			void build_test_data_lod1(Builder &builder) {

				const size_t expected_num_vertices = 16;
				const size_t expected_num_faces    = 5;


				// Start surface building.
				m_index_counter = 0;
				builder.begin_surface(expected_num_vertices, expected_num_faces);


				// Add triangle.
				add_triangle(Point( 0.0, 0.0, 0.0), Point( 1.0, 0.0, 0.0), Point( 0.0,  1.0, 0.0), generate_random_color(), builder);

		        // Add the same triangle.
				add_triangle(Point( 0.0, 0.0, 0.0), Point( 1.0, 0.0, 0.0), Point( 0.0,  1.0, 0.0), generate_random_color(), builder);		        

		        // Add triangle with one common edge.
		        add_triangle(Point(-1.0, 0.0, 0.0), Point( 0.0, 0.0, 0.0), Point( 0.0,  0.5, 0.0), generate_random_color(), builder);		        

		        // Add triangle with one common vertex.
		        add_triangle(Point( 0.0, 0.0, 0.0), Point(-1.0, 0.0, 0.0), Point( 0.0, -1.0, 0.0), generate_random_color(), builder);	

		        // Add quad.
		        add_quad(Point(1.0, 0.0, 0.0), Point(2.0, 0.0, 0.0), Point(2.0, 1.0, 0.0), Point(1.0, 1.0, 0.0), generate_random_color(), builder);


		        // End surface building.
		        builder.end_surface();
			}

			void add_triangle(const Point &a, const Point &b, const Point &c, const Color &color, Builder &builder) {

				// Add vertices of the triangle facet.
		        builder.add_vertex(a);
		        builder.add_vertex(b);
		        builder.add_vertex(c);

		        // Add new triangle facet.
		        const Color_facet_handle cfh = builder.begin_facet();

		        builder.add_vertex_to_facet(m_index_counter++);
		        builder.add_vertex_to_facet(m_index_counter++);
		        builder.add_vertex_to_facet(m_index_counter++);
		        
		        builder.end_facet();
		        
		        // Add triangle facet color.
		        m_facet_colors[cfh] = color;
			}

			void add_quad(const Point &a, const Point &b, const Point &c, const Point &d, const Color &color, Builder &builder) {

				// Add vertices of the quad facet.
		        builder.add_vertex(a);
		        builder.add_vertex(b);
		        builder.add_vertex(c);
		        builder.add_vertex(d);

		        // Add facets.
		        const Color_facet_handle cfh = builder.begin_facet();

		        builder.add_vertex_to_facet(m_index_counter++);
		        builder.add_vertex_to_facet(m_index_counter++);
		        builder.add_vertex_to_facet(m_index_counter++);
		        builder.add_vertex_to_facet(m_index_counter++);

		        builder.end_facet();

		        // Add quad facet color.
		        m_facet_colors[cfh] = color;
			}

			Color generate_random_color() {

				const FT r = m_rand.get_int(0, 256);
				const FT g = m_rand.get_int(0, 256);
				const FT b = m_rand.get_int(0, 256);

				return Color(r, g, b);
			}

			void build_from_cdt_and_buildings_lod1(Builder &builder) {

				const size_t expected_num_vertices = estimate_number_of_vertices_lod1();
				const size_t expected_num_faces    = estimate_number_of_faces_lod1();


				// Start surface building.
				m_index_counter = 0;
				builder.begin_surface(expected_num_vertices, expected_num_faces);


				// Add all buildings.
				for (Building_iterator bit = m_buildings.begin(); bit != m_buildings.end(); ++bit) {
				
					const auto &building = (*bit).second;
					add_new_building_lod1(building, builder);
				}


				// Add ground.
				if (!m_ground_set) estimate_ground();
				add_ground(builder);


				// End surface building.
				builder.end_surface();
			}

			// Improve this function.
			size_t estimate_number_of_vertices_lod0() {
				return m_buildings.size() * 4 + 4;
			}

			// Improve this function.
			size_t estimate_number_of_faces_lod0() {
				return m_buildings.size() + 1;
			}

			// Improve this function.
			size_t estimate_number_of_vertices_lod1() {
				return m_buildings.size() * 4 * 2 + 4;
			}

			// Improve this function.
			size_t estimate_number_of_faces_lod1() {
				return m_buildings.size() * 6 + 1;
			}

			template<class BuildingTmp>
			void add_new_building_lod0(const BuildingTmp &building, Builder &builder) {
				
				const auto &boundaries = building.boundaries;
				const size_t num_boundaries = boundaries.size();

				if (num_boundaries == 1) {
					
					add_building_structure_from_one_boundary_lod0(building, builder);
					return;
				}

				if (num_boundaries > 1) {

					add_building_structure_from_multiple_boundaries_lod0(building, builder);
					return;
				}

				if (num_boundaries == 0 && !building.faces.empty() && !m_use_boundaries) {

					add_building_structure_from_one_boundary_lod0(building, builder);
					return;
				}

				if (num_boundaries == 0) return;
			}

			template<class BuildingTmp>
			void add_new_building_lod1(const BuildingTmp &building, Builder &builder) {
				
				const auto &boundaries = building.boundaries;
				const size_t num_boundaries = boundaries.size();

				if (num_boundaries == 1) {
					
					add_building_structure_from_one_boundary_lod1(building, builder);
					return;
				}

				if (num_boundaries > 1) {

					add_building_structure_from_multiple_boundaries_lod1(building, builder);
					return;
				}

				if (num_boundaries == 0 && !building.faces.empty() && !m_use_boundaries) {

					add_building_structure_from_one_boundary_lod1(building, builder);
					return;
				}

				if (num_boundaries == 0) return;
			}

			template<class BuildingTmp>
			void add_building_structure_from_one_boundary_lod0(const BuildingTmp &building, Builder &builder) {

				const Color color = building.color;
				if (m_use_boundaries) {

					const auto &boundary = building.boundaries[0];
					add_horizontal_boundary(boundary, color, FT(0), builder); // floor

				} else {

					const auto &faces = building.faces;
					add_horizontal_triangulation(faces, color, FT(0), builder); // floor
				}
			}

			template<class BuildingTmp>
			void add_building_structure_from_one_boundary_lod1(const BuildingTmp &building, Builder &builder) {

				const FT height   = building.height;
				const Color color = building.color;

				if (m_use_boundaries) {
	
					const auto &boundary = building.boundaries[0];

					add_horizontal_boundary(boundary, color, FT(0) , builder); // floor
					add_horizontal_boundary(boundary, color, height, builder); // roof

					add_walls_from_boundary(boundary, color, FT(0), height, builder); // walls

				} else {

					const auto &faces = building.faces;

					add_horizontal_triangulation(faces, color, FT(0) , builder); // floor
					add_horizontal_triangulation(faces, color, height, builder); // roof

					add_walls_from_triangles(faces, color, FT(0), height, builder); // walls
				}
			}

			template<class BoundaryTmp>
			void add_horizontal_boundary(const BoundaryTmp &boundary, const Color &color, const FT height, Builder &builder) {

				Point point;

				// Add vertices.
				const size_t num_vertices = boundary.size();
				for (size_t i = 0; i < num_vertices; ++i) {

					const auto &vertex_handle = boundary[i];
					const auto &vertex = vertex_handle->point();

					point = Point(vertex.x(), vertex.y(), height);
					builder.add_vertex(point);
				}

				// Add new polygonal facet.
				const Color_facet_handle cfh = builder.begin_facet();
				for (size_t i = 0; i < num_vertices; ++i) builder.add_vertex_to_facet(m_index_counter++);
				builder.end_facet();

				// Add polygonal facet color.
		        m_facet_colors[cfh] = color;
			}

			template<class BuildingTmp>
			void add_building_structure_from_multiple_boundaries_lod0(const BuildingTmp &building, Builder &builder) {

				const Color color = building.color;
				const auto &faces = building.faces;

				add_horizontal_triangulation(faces, color, FT(0), builder); // floor
			}

			template<class BuildingTmp>
			void add_building_structure_from_multiple_boundaries_lod1(const BuildingTmp &building, Builder &builder) {

				const Color color = building.color;
				const FT height   = building.height;

				const auto &faces = building.faces;

				add_horizontal_triangulation(faces, color, FT(0) , builder); // floor
				add_horizontal_triangulation(faces, color, height, builder); // roof

				// Add walls.
				if (m_use_boundaries) {

					const size_t num_boundaries = building.boundaries.size();
					for (size_t i = 0; i < num_boundaries; ++i) {

						const auto &boundary = building.boundaries[i];
						add_walls_from_boundary(boundary, color, FT(0), height, builder);		
					}

				} else add_walls_from_triangles(faces, color, FT(0), height, builder);
			}

			template<class FaceHandlesTmp>
			void add_horizontal_triangulation(const FaceHandlesTmp &fhs, const Color &color, const FT height, Builder &builder) {

				// Add all faces.
				const size_t num_face_handles = fhs.size();
				for (size_t i = 0; i < num_face_handles; ++i) {

					const auto &fh = fhs[i];
					const auto &triangle = m_cdt.triangle(fh);

					add_triangle_face(triangle, color, height, builder);
				}
			}

			template<class TriangleFaceTmp>
			void add_triangle_face(const TriangleFaceTmp &triangle, const Color &color, const FT height, Builder &builder) {

				const auto &va = triangle.vertex(0);
				const auto &vb = triangle.vertex(1);
				const auto &vc = triangle.vertex(2);

				const Point a = Point(va.x(), va.y(), height);
				const Point b = Point(vb.x(), vb.y(), height);
				const Point c = Point(vc.x(), vc.y(), height);

				add_triangle(a, b, c, color, builder);
			}

			template<class BoundaryTmp>
			void add_walls_from_boundary(const BoundaryTmp &boundary, const Color &color, const FT height_floor, const FT height_roof, Builder &builder) {

				const size_t num_vertices = boundary.size();
				for (size_t i = 0; i < num_vertices; ++i) {

					const size_t ip = (i + 1) % num_vertices;
					add_quadrilateral_wall(boundary[i], boundary[ip], color, height_floor, height_roof, builder);
				}
			}

			template<class FaceHandlesTmp>
			void add_walls_from_triangles(const FaceHandlesTmp &fhs, const Color &color, const FT height_floor, const FT height_roof, Builder &builder) {

				const size_t num_face_handles = fhs.size();
				for (size_t i = 0; i < num_face_handles; ++i) {
					
					const auto &fh = fhs[i];
					for (size_t j = 0; j < 3; ++j) {

						const size_t jp = (j + 1) % 3;
						add_quadrilateral_wall(fh->vertex(j), fh->vertex(jp), color, height_floor, height_roof, builder);		
					}
				}
			}

			template<class VertexHandleTmp>
			void add_quadrilateral_wall(const VertexHandleTmp &vi, const VertexHandleTmp &vj, const Color &color, const FT height_floor, const FT height_roof, Builder &builder) {

				const auto &va = vi->point();
				const auto &vb = vj->point();

				const Point a = Point(va.x(), va.y(), height_floor);
				const Point b = Point(vb.x(), vb.y(), height_floor);

				const Point c = Point(vb.x(), vb.y(), height_roof);
				const Point d = Point(va.x(), va.y(), height_roof);

				add_quad(a, b, c, d, color, builder);
			}

			void estimate_ground() {

				const FT big_value = FT(1000000000);
				FT min_x = big_value, min_y = big_value, max_x = -big_value, max_y = -big_value;

				for (auto vit = m_cdt.finite_vertices_begin(); vit != m_cdt.finite_vertices_end(); ++vit) {
					const auto &p = vit->point();

					min_x = CGAL::min(min_x, p.x());
					min_y = CGAL::min(min_y, p.y());

					max_x = CGAL::max(max_x, p.x());
					max_y = CGAL::max(max_y, p.y());
				}
				
				m_ground.clear();
				m_ground.resize(4);

				m_ground[0] = Point(min_x, min_y, FT(0));
				m_ground[1] = Point(max_x, min_y, FT(0));
				m_ground[2] = Point(max_x, max_y, FT(0));
				m_ground[3] = Point(min_x, max_y, FT(0));
			}

			void add_ground(Builder &builder) {

				assert(!m_ground.empty());

				const size_t num_vertices = m_ground.size();
				assert(num_vertices == 4);

				const Point &a = m_ground[3];
				const Point &b = m_ground[2];
				const Point &c = m_ground[1];
				const Point &d = m_ground[0];

				const Color color(169, 169, 169);
				add_quad(a, b, c, d, color, builder);
			}
		}; 


		// The main LOD1 reconstruction class.
		template<class KernelTraits, class CDTInput, class BuildingsInput, class MeshOutput>
		class Level_of_detail_reconstruction {

		public:
			typedef KernelTraits   Kernel;
			typedef CDTInput 	   CDT;
			typedef BuildingsInput Buildings;
			typedef MeshOutput 	   Mesh;

			// Extra.
			using HDS = typename Mesh::HalfedgeDS;

			using Mesh_facet_handle = typename Mesh::Facet_const_handle;
			using Mesh_builder 		= Build_mesh<Kernel, HDS, CDT, Buildings, Mesh_facet_handle>;
			using Mesh_facet_colors = typename Mesh_builder::Facet_colors;

			using Point  = typename Mesh_builder::Point;
			using Ground = typename Mesh_builder::Ground;

			Level_of_detail_reconstruction() : m_use_boundaries(true) { }

			void use_boundaries(const bool new_state) {
				m_use_boundaries = new_state;
			}

			void reconstruct_lod0(const CDT &cdt, const Buildings &buildings, const Ground &ground, Mesh &mesh, Mesh_facet_colors &mesh_facet_colors) const {

				Mesh_builder mesh_builder(cdt, buildings, mesh_facet_colors);
				
				mesh_builder.set_ground(ground);
				mesh_builder.use_boundaries(m_use_boundaries);

				mesh_builder.set_builder_type(Mesh_builder::Builder_type::LOD0);
				mesh.delegate(mesh_builder);
			}

			void reconstruct_lod1(const CDT &cdt, const Buildings &buildings, const Ground &ground, Mesh &mesh, Mesh_facet_colors &mesh_facet_colors) const {
			
				Mesh_builder mesh_builder(cdt, buildings, mesh_facet_colors);

				mesh_builder.set_ground(ground);
				mesh_builder.use_boundaries(m_use_boundaries);

				mesh_builder.set_builder_type(Mesh_builder::Builder_type::LOD1);
				mesh.delegate(mesh_builder);
			}

		private:
			bool m_use_boundaries;
		};
	}
}

#endif // CGAL_LEVEL_OF_DETAIL_RECONSTRUCTION_H