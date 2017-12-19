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
#include <CGAL/Utils/Level_of_detail_utils.h>

namespace CGAL {

	namespace LOD {

		// Mesh builder for the reconstruction below.
		template<class KernelTraits, class HDSInput, class CDTInput, class BuildingsInput, class FacetHandle>
		class Build_mesh : public Modifier_base<HDSInput> {
		
		public:
			typedef KernelTraits 	 Kernel;
    		typedef HDSInput 		 HDS;
    		typedef CDTInput 		 CDT;
    		typedef BuildingsInput 	 Buildings;
    		typedef FacetHandle 	 Color_facet_handle;
    		typedef FacetHandle 	 Facet_handle;

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

    		typedef Level_of_detail_utils_simple<Kernel> Simple_utils;
			typename Kernel::Compute_squared_distance_2 squared_distance;

    		enum class Builder_type { LOD0, LOD1, ROOFS, WALLS };


    		Build_mesh(const CDT &cdt, const Buildings &buildings, Facet_colors &facet_colors) : 
    		m_build_type(Build_type::CDT_AND_BUILDINGS), 
    		m_cdt(cdt), 
    		m_buildings(buildings), 
    		m_facet_colors(facet_colors), 
    		m_index_counter(0), 
    		m_ground_set(false), 
    		m_builder_type(Builder_type::LOD1),
    		m_use_boundaries(true),
    		m_height_threshold(FT(1) / FT(1000000)),
    		m_area_threshold(FT(1) / FT(1000)),
    		m_num_roofs(-1),
    		m_num_walls(-1) { }

			void operator()(HDS &hds) {

				Builder builder(hds, false);
				switch (m_builder_type) {

					case Builder_type::LOD0:
						build_lod0(builder);
						break;

					case Builder_type::LOD1:
						build_lod1(builder);
						break;

					case Builder_type::ROOFS:
						build_lod1_roofs(builder);
						break;

					case Builder_type::WALLS:
						build_lod1_walls(builder);
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


			int get_number_of_roofs() const {

				assert(m_num_roofs >= 0);
				return m_num_roofs;
			}

			int get_number_of_walls() const {

				assert(m_num_walls >= 0);
				return m_num_walls;
			}

			void get_roofs_faces(std::vector< std::vector<Facet_handle> > &roofs_faces) const {
				roofs_faces = m_roofs_faces;
			}

			void get_walls_faces(std::vector< std::vector<Facet_handle> > &walls_faces) const {
				walls_faces = m_walls_faces;
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

			const FT 	 m_height_threshold;
			const FT 	 m_area_threshold;
			Simple_utils m_simple_utils;

			int m_num_roofs;
			int m_num_walls;

			std::vector< std::vector<Facet_handle> > m_roofs_faces, m_walls_faces;


			// EXTRA CODE!

			void build_lod1_roofs(Builder &builder) {

				const size_t expected_num_vertices = 10;
				const size_t expected_num_faces    = 10;

				size_t num_roofs 	 = 0;
				size_t index_counter = 0;

				const size_t num_buildings = m_buildings.size();

				m_roofs_faces.clear();
				m_roofs_faces.resize(num_buildings);

				size_t building_index = 0;
				builder.begin_surface(expected_num_vertices, expected_num_faces);
				for (Building_iterator bit = m_buildings.begin(); bit != m_buildings.end(); ++bit, ++building_index) {
				
					const auto &building = (*bit).second;
					if (!is_valid_building(building)) continue;

					if (building.is_oriented) {

						std::cerr << std::endl << "ERROR: Oriented method is not working! Please turn off the quality option!" << std::endl << std::endl;
						exit(EXIT_FAILURE);
					}
					const FT height = building.height;
					++num_roofs;
	
					const auto &faces = building.faces;
					add_horizontal_triangulation_custom(faces, height, builder, index_counter, m_roofs_faces[building_index]);
				}
				
				builder.end_surface();
				assert(m_num_roofs == static_cast<int>(num_roofs));
			}

			template<class FaceHandlesTmp>
			void add_horizontal_triangulation_custom(
				const FaceHandlesTmp &fhs, const FT height, 
				Builder &builder, size_t &index_counter, std::vector<Facet_handle> &faces) {

				const size_t num_face_handles = fhs.size();
				for (size_t i = 0; i < num_face_handles; ++i) {

					const auto &fh = fhs[i];
					const auto &triangle = m_cdt.triangle(fh);

					add_triangle_face_custom(triangle, height, builder, index_counter, faces);
				}
			}

			template<class TriangleFaceTmp>
			void add_triangle_face_custom(
				const TriangleFaceTmp &triangle, const FT height, 
				Builder &builder, size_t &index_counter, std::vector<Facet_handle> &faces) {

				const auto &va = triangle.vertex(0);
				const auto &vb = triangle.vertex(1);
				const auto &vc = triangle.vertex(2);

				const Point a = Point(va.x(), va.y(), height);
				const Point b = Point(vb.x(), vb.y(), height);
				const Point c = Point(vc.x(), vc.y(), height);

				add_triangle_custom(a, b, c, builder, index_counter, faces);
			}

			void add_triangle_custom(
				const Point &a, const Point &b, const Point &c, 
				Builder &builder, size_t &index_counter, std::vector<Facet_handle> &faces) {

		        builder.add_vertex(a);
		        builder.add_vertex(b);
		        builder.add_vertex(c);

		        const Facet_handle fh = builder.begin_facet();

		        builder.add_vertex_to_facet(index_counter++);
		        builder.add_vertex_to_facet(index_counter++);
		        builder.add_vertex_to_facet(index_counter++);
		        
		        builder.end_facet();
		      	faces.push_back(fh);
			}

			void build_lod1_walls(Builder &builder) {

				const size_t expected_num_vertices = 10;
				const size_t expected_num_faces    = 10;

				size_t num_walls 	 = 0;
				size_t index_counter = 0;

				const size_t num_buildings = m_buildings.size();
				
				m_walls_faces.clear();
				m_walls_faces.resize(num_buildings);

				size_t building_index = 0;
				builder.begin_surface(expected_num_vertices, expected_num_faces);
				for (Building_iterator bit = m_buildings.begin(); bit != m_buildings.end(); ++bit, ++building_index) {
				
					const auto &building = (*bit).second;
					if (!is_valid_building(building)) continue;

					if (building.is_oriented) {

						std::cerr << std::endl << "ERROR: Oriented method is not working! Please turn off the quality option!" << std::endl << std::endl;
						exit(EXIT_FAILURE);
					}
					
					const FT height = building.height;
					assert(building.boundaries.size() > 0);

					const auto &boundary = building.boundaries[0];
					add_walls_from_unoriented_boundary_custom(boundary, FT(0), height, builder, num_walls, index_counter, m_walls_faces[building_index]);
				}

				builder.end_surface();
				assert(m_num_walls == static_cast<int>(num_walls));
			}

			template<class BoundaryTmp>
			void add_walls_from_unoriented_boundary_custom(
				const BoundaryTmp &boundary, const FT height_floor, const FT height_roof, 
				Builder &builder, size_t &num_walls, size_t &index_counter, std::vector<Facet_handle> &faces) {
				
				const size_t num_vertices = boundary.size();
				for (size_t i = 0; i < num_vertices; i += 2) {
					
					const size_t ip = i + 1;
					assert(ip < num_vertices);

					++num_walls;
					add_quadrilateral_wall_custom(boundary[i], boundary[ip], height_floor, height_roof, builder, index_counter, faces);
				}
			}

			template<class VertexHandleTmp>
			void add_quadrilateral_wall_custom(
				const VertexHandleTmp &vi, const VertexHandleTmp &vj, const FT height_floor, const FT height_roof, 
				Builder &builder, size_t &index_counter, std::vector<Facet_handle> &faces) {

				const auto &va = vi->point();
				const auto &vb = vj->point();

				const Point a = Point(va.x(), va.y(), height_floor);
				const Point b = Point(vb.x(), vb.y(), height_floor);

				const Point c = Point(vb.x(), vb.y(), height_roof);
				const Point d = Point(va.x(), va.y(), height_roof);

				add_quad_custom(a, b, c, d, builder, index_counter, faces);
			}

			void add_quad_custom(
				const Point &a, const Point &b, const Point &c, const Point &d,
				Builder &builder, size_t &index_counter, std::vector<Facet_handle> &faces) {

		        builder.add_vertex(a);
		        builder.add_vertex(b);
		        builder.add_vertex(c);
		        builder.add_vertex(d);

		        const Facet_handle fh = builder.begin_facet();

		        builder.add_vertex_to_facet(index_counter++);
		        builder.add_vertex_to_facet(index_counter++);
		        builder.add_vertex_to_facet(index_counter++);
		        builder.add_vertex_to_facet(index_counter++);

		        builder.end_facet();
		        faces.push_back(fh);
			}






			// MAIN CODE!

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
					if (!is_valid_building(building)) continue;

					if (building.is_oriented) add_new_building_lod0(building, builder); // default is oriented
					else add_new_building_lod0_unoriented(building, builder);
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

				return Color(CGAL::to_double(r), CGAL::to_double(g), CGAL::to_double(b));
			}

			void build_from_cdt_and_buildings_lod1(Builder &builder) {

				const size_t expected_num_vertices = estimate_number_of_vertices_lod1();
				const size_t expected_num_faces    = estimate_number_of_faces_lod1();

				// Start surface building.
				m_index_counter = 0;
				builder.begin_surface(expected_num_vertices, expected_num_faces);

				// Add all buildings.
				m_num_roofs = 0;
				m_num_walls = 0;

				for (Building_iterator bit = m_buildings.begin(); bit != m_buildings.end(); ++bit) {
				
					const auto &building = (*bit).second;
					if (!is_valid_building(building)) continue;

					++m_num_roofs;

					if (building.is_oriented) {
						add_new_building_lod1(building, builder); // default is oriented
						std::cerr << std::endl << "WARNING: number of walls will be wrong! See complexity." << std::endl << std::endl;
					}
					else add_new_building_lod1_unoriented(building, builder);
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
			bool is_valid_building(const BuildingTmp &building) {

				const auto &faces = building.faces;
				if (faces.size() < 2) return false;

				const FT height = building.height;
				if (height < m_height_threshold) return false;

				return true;
			}

			FT compute_ground_area() {

				const Point &a = m_ground[0];
				const Point &b = m_ground[1];
				const Point &c = m_ground[2];

				const FT width  = static_cast<FT>(CGAL::sqrt(CGAL::to_double(squared_distance(a, b))));
				const FT height = static_cast<FT>(CGAL::sqrt(CGAL::to_double(squared_distance(b, c))));

				return width * height;
			}

			template<class BuildingTmp>
			FT compute_building_area(const BuildingTmp &building) {

				const auto &faces = building.faces;

				FT building_area = FT(0);
				for (size_t i = 0; i < faces.size(); ++i) {
				
					const auto &triangle = m_cdt.triangle(faces[i]);
					building_area += triangle.area();
				}

				return building_area;
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
			void add_new_building_lod0_unoriented(const BuildingTmp &building, Builder &builder) {
				
				const auto &boundaries = building.boundaries;
				const size_t num_boundaries = boundaries.size();

				assert(num_boundaries > 0);
				if (num_boundaries <= 0) {

					std::cerr << "Error: num_boundaries <= 0, add_new_building_lod0_unoriented function reconstruction!" << std::endl;
					exit(EXIT_FAILURE);
				}

				add_unoriented_building_structure_lod0(building, builder);
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
			void add_new_building_lod1_unoriented(const BuildingTmp &building, Builder &builder) {

				const auto &boundaries = building.boundaries;
				const size_t num_boundaries = boundaries.size();

				assert(num_boundaries > 0);
				if (num_boundaries <= 0) {

					std::cerr << "Error: num_boundaries <= 0, add_new_building_lod1_unoriented function reconstruction!" << std::endl;
					exit(EXIT_FAILURE);
				}

				add_unoriented_building_structure_lod1(building, builder);
			}

			template<class BuildingTmp>
			void add_building_structure_from_one_boundary_lod0(const BuildingTmp &building, Builder &builder) {

				const Color color = building.color;
				if (m_use_boundaries) {

					const auto &boundary = building.boundaries[0];
					add_horizontal_boundary(boundary, color, FT(1) / FT(1000), builder); // floor

				} else {

					const auto &faces = building.faces;
					add_horizontal_triangulation(faces, color, FT(1) / FT(1000), builder); // floor
				}
			}

			template<class BuildingTmp>
			void add_unoriented_building_structure_lod0(const BuildingTmp &building, Builder &builder) {

				const Color color = building.color;
				const auto &faces = building.faces;

				add_horizontal_triangulation(faces, color, FT(1) / FT(1000), builder); // floor
			}

			template<class BuildingTmp>
			void add_building_structure_from_one_boundary_lod1(const BuildingTmp &building, Builder &builder) {

				const FT height   = building.height;
				const Color color = building.color;

				if (m_use_boundaries) {
	
					const auto &boundary = building.boundaries[0];

					add_horizontal_boundary(boundary, color, height, builder); 		  // roof
					add_walls_from_boundary(boundary, color, FT(0), height, builder); // walls

				} else {

					const auto &faces = building.faces;

					add_horizontal_triangulation(faces, color, height, builder);  	// roof
					add_walls_from_triangles(faces, color, FT(0), height, builder); // walls
				}
			}

			template<class BuildingTmp>
			void add_unoriented_building_structure_lod1(const BuildingTmp &building, Builder &builder) {

				const FT height   = building.height;
				const Color color = building.color;

				const auto &faces = building.faces;
				add_horizontal_triangulation(faces, color, height, builder); // roof
	
				const auto &boundary = building.boundaries[0];
				add_walls_from_unoriented_boundary(boundary, color, FT(0), height, builder); // walls
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

				add_horizontal_triangulation(faces, color, FT(1) / FT(1000), builder); // floor
			}

			template<class BuildingTmp>
			void add_building_structure_from_multiple_boundaries_lod1(const BuildingTmp &building, Builder &builder) {

				const Color color = building.color;
				const FT height   = building.height;

				const auto &faces = building.faces;

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

			template<class BoundaryTmp>
			void add_walls_from_unoriented_boundary(const BoundaryTmp &boundary, const Color &color, const FT height_floor, const FT height_roof, Builder &builder) {
				
				const size_t num_vertices = boundary.size();
				for (size_t i = 0; i < num_vertices; i += 2) {
					
					const size_t ip = i + 1;
					assert(ip < num_vertices);

					++m_num_walls;
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

				if (num_vertices != 4) {
					std::cerr << "Error: num_vertices != 4, add_ground function reconstruction!" << std::endl;
					exit(EXIT_FAILURE);
				}

				const Point &a = m_ground[3];
				const Point &b = m_ground[2];
				const Point &c = m_ground[1];
				const Point &d = m_ground[0];

				const Color color(169, 169, 169);
				add_quad(a, b, c, d, color, builder);
			}
		}; 


		// The main reconstruction class.
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

			Level_of_detail_reconstruction() : 
			m_use_boundaries(true),
			m_num_roofs(-1),
			m_num_walls(-1) { }

			void use_boundaries(const bool new_state) {
				m_use_boundaries = new_state;
			}

			int get_number_of_roofs() const {
				
				assert(m_num_roofs >= 0);
				return m_num_roofs;
			}

			int get_number_of_walls() const {
				
				assert(m_num_walls >= 0);
				return m_num_walls;
			}

			void reconstruct_lod0(const CDT &cdt, const Buildings &buildings, const Ground &ground, Mesh &mesh, Mesh_facet_colors &mesh_facet_colors) {

				Mesh_builder mesh_builder(cdt, buildings, mesh_facet_colors);
				
				mesh_builder.set_ground(ground);
				mesh_builder.use_boundaries(m_use_boundaries);

				mesh_builder.set_builder_type(Mesh_builder::Builder_type::LOD0);
				mesh.delegate(mesh_builder);
			}

			void reconstruct_lod1(const CDT &cdt, const Buildings &buildings, const Ground &ground, Mesh &mesh, Mesh_facet_colors &mesh_facet_colors) {

				Mesh_builder mesh_builder(cdt, buildings, mesh_facet_colors);

				mesh_builder.set_ground(ground);
				mesh_builder.use_boundaries(m_use_boundaries);

				mesh_builder.set_builder_type(Mesh_builder::Builder_type::LOD1);
				mesh.delegate(mesh_builder);

				set_lod1_metrics(mesh_builder);

				reconstruct_lod1_roofs(cdt, buildings);
				reconstruct_lod1_walls(cdt, buildings);
			}

			void get_roofs(Mesh &roofs) const {
				roofs = m_roofs;
			}

			void get_walls(Mesh &walls) const {
				walls = m_walls;
			}

			void get_roofs_faces(std::vector< std::vector<Mesh_facet_handle> > &roofs_faces) const {
				roofs_faces = m_roofs_faces;
			}

			void get_walls_faces(std::vector< std::vector<Mesh_facet_handle> > &walls_faces) const {
				walls_faces = m_walls_faces;
			}

		private:

			bool m_use_boundaries;

			int m_num_roofs;
			int m_num_walls;

			Mesh m_roofs;
			Mesh m_walls;

			std::vector< std::vector<Mesh_facet_handle> > m_roofs_faces, m_walls_faces;


			void reconstruct_lod1_roofs(const CDT &cdt, const Buildings &buildings) {
				set_lod1_roofs(cdt, buildings);
			}

			void reconstruct_lod1_walls(const CDT &cdt, const Buildings &buildings) {
				set_lod1_walls(cdt, buildings);
			}

			void set_lod1_metrics(const Mesh_builder &mesh_builder) {

				m_num_roofs = mesh_builder.get_number_of_roofs();
				m_num_walls = mesh_builder.get_number_of_walls();
			}

			void set_lod1_roofs(const CDT &cdt, const Buildings &buildings) {

				Mesh_facet_colors stub;
				Mesh_builder roofs_builder(cdt, buildings, stub);

				roofs_builder.set_builder_type(Mesh_builder::Builder_type::ROOFS);

				m_roofs.clear();
				m_roofs.delegate(roofs_builder);

				roofs_builder.get_roofs_faces(m_roofs_faces);
			}

			void set_lod1_walls(const CDT &cdt, const Buildings &buildings) {

				Mesh_facet_colors stub;
				Mesh_builder walls_builder(cdt, buildings, stub);

				walls_builder.set_builder_type(Mesh_builder::Builder_type::WALLS);

				m_walls.clear();
				m_walls.delegate(walls_builder);

				walls_builder.get_walls_faces(m_walls_faces);
			}
		};
	}
}

#endif // CGAL_LEVEL_OF_DETAIL_RECONSTRUCTION_H