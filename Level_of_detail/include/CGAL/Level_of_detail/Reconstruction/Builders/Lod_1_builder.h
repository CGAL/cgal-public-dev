#ifndef CGAL_LEVEL_OF_DETAIL_LOD_1_BUILDER_H
#define CGAL_LEVEL_OF_DETAIL_LOD_1_BUILDER_H

// STL includes.
#include <vector>

// CGAL includes.
#include <CGAL/Modifier_base.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>

namespace CGAL {

	namespace Level_of_detail {

        template<
        class InputKernel, 
        class InputHDS, 
        class InputBuilding, 
        class InputBuildings, 
        class InputGround,
        class GroundFace, 
        class WallFaces,
        class RoofFaces>
		class Lod_1_builder : public Modifier_base<InputHDS> {
		
		public:
            using Kernel   = InputKernel;
            using HDS      = InputHDS;
            using Building = InputBuilding;

            using Buildings = InputBuildings;
            using Ground    = InputGround;

            using Ground_face = GroundFace;
            using Wall_faces  = WallFaces;
            using Roof_faces  = RoofFaces;

            using FT      = typename Kernel::FT;
            using Point_2 = typename Kernel::Point_2;
            using Point_3 = typename Kernel::Point_3;

            using Floor_edges = typename Building::Floor_edges;
            using Floor_faces = typename Building::Floor_faces;

            using Ground_points_iterator     = typename Ground::const_iterator;
            using Const_floor_edges_iterator = typename Building::Const_floor_edges_iterator;
            using Const_floor_faces_iterator = typename Building::Const_floor_faces_iterator;
            using Const_buildings_iterator   = typename Buildings::const_iterator;

            using Points_3 = std::vector<Point_3>;
            using Builder  = CGAL::Polyhedron_incremental_builder_3<HDS>;
            
            Lod_1_builder(const Buildings &buildings, const Ground &ground, Ground_face &ground_face, Wall_faces &wall_faces, Roof_faces &roof_faces) :
            m_buildings(buildings),
            m_ground(ground),
            m_ground_face(ground_face),
            m_wall_faces(wall_faces),
            m_roof_faces(roof_faces),
            m_index_counter(0)
            { }

			void operator()(HDS &hds) {
				
				Builder builder(hds, false);

				const size_t expected_num_vertices = estimate_number_of_vertices();
				const size_t expected_num_faces    = estimate_number_of_faces();

				m_index_counter = 0;
				builder.begin_surface(expected_num_vertices, expected_num_faces);

                m_wall_faces.clear();
                m_roof_faces.clear();

				for (Const_buildings_iterator cb_it = m_buildings.begin(); cb_it != m_buildings.end(); ++cb_it) {
				    const Building &building = *cb_it;
				    
                    if (!building.is_valid()) continue;
				    add_new_building(building, builder);
				}

				add_ground(builder);
				builder.end_surface();
			}

        private:
            const Buildings &m_buildings;
            const Ground    &m_ground;
            
            Ground_face &m_ground_face;
            Wall_faces  &m_wall_faces;
            Roof_faces  &m_roof_faces;

            size_t m_index_counter;

            inline size_t estimate_number_of_vertices() {
				return m_buildings.size() * 8 + 4;
			}

			inline size_t estimate_number_of_faces() {
				return m_buildings.size() * 6 + 1;
			}

			void add_new_building(const Building &building, Builder &builder) {

                const FT building_height       = building.height();
				const Floor_edges &floor_edges = building.floor_edges();
                const Floor_faces &floor_faces = building.floor_faces();
                
				add_wall_faces(floor_edges, building_height, builder);
                add_roof_faces(floor_faces, building_height, builder);
			}

            void add_wall_faces(const Floor_edges &floor_edges, const FT building_height, Builder &builder) {

                size_t i = 0;
                const FT base_z = m_ground.begin()->z();

                Points_3 face_vertices(4);
                for (Const_floor_edges_iterator segment = floor_edges.begin(); segment != floor_edges.end(); ++segment, ++i) {
                    
                    const Point_2 &p1 = segment->source();
                    const Point_2 &p2 = segment->target();

                    const Point_3 a = Point_3(p1.x(), p1.y(), base_z);
                    const Point_3 b = Point_3(p2.x(), p2.y(), base_z);
                    const Point_3 c = Point_3(p2.x(), p2.y(), base_z + building_height);
                    const Point_3 d = Point_3(p1.x(), p1.y(), base_z + building_height);

                    face_vertices[0] = a;
                    face_vertices[1] = b;
                    face_vertices[2] = c;
                    face_vertices[3] = d;
                    m_wall_faces.push_back(face_vertices);

                    add_quadrilateral(a, b, c, d, builder);
                }
            }

            void add_roof_faces(const Floor_faces &floor_faces, const FT building_height, Builder &builder) {

                size_t i = 0;
                const FT base_z = m_ground.begin()->z();

                Points_3 face_vertices(3);
                for (Const_floor_faces_iterator triangle = floor_faces.begin(); triangle != floor_faces.end(); ++triangle, ++i) {
                    
                    const Point_2 &p1 = triangle->vertex(0);
                    const Point_2 &p2 = triangle->vertex(1);
                    const Point_2 &p3 = triangle->vertex(2);

                    const Point_3 a = Point_3(p1.x(), p1.y(), base_z + building_height);
                    const Point_3 b = Point_3(p2.x(), p2.y(), base_z + building_height);
                    const Point_3 c = Point_3(p3.x(), p3.y(), base_z + building_height);

                    face_vertices[0] = a;
                    face_vertices[1] = b;
                    face_vertices[2] = c;
                    m_roof_faces.push_back(face_vertices);

                    add_triangle(a, b, c, builder);
                }
            }

            void add_triangle(const Point_3 &a, const Point_3 &b, const Point_3 &c, Builder &builder) {

		        builder.add_vertex(a);
		        builder.add_vertex(b);
		        builder.add_vertex(c);

		        builder.begin_facet();
		        builder.add_vertex_to_facet(m_index_counter++);
		        builder.add_vertex_to_facet(m_index_counter++);
		        builder.add_vertex_to_facet(m_index_counter++);
		        builder.end_facet();
			}

            void add_ground(Builder &builder) {
				
				const size_t num_vertices = m_ground.size();
                CGAL_precondition(num_vertices == 4);

                m_ground_face.clear();
                m_ground_face.resize(num_vertices);

                size_t i = 0;
                for (Ground_points_iterator point = m_ground.begin(); point != m_ground.end(); ++point, ++i)
                    m_ground_face[i] = *point;

				const Point_3 &a = m_ground_face[0];
				const Point_3 &b = m_ground_face[1];
				const Point_3 &c = m_ground_face[2];
				const Point_3 &d = m_ground_face[3];

				add_quadrilateral(a, b, c, d, builder);
			}

            void add_quadrilateral(const Point_3 &a, const Point_3 &b, const Point_3 &c, const Point_3 &d, Builder &builder) {

		        builder.add_vertex(a);
		        builder.add_vertex(b);
		        builder.add_vertex(c);
		        builder.add_vertex(d);

		        builder.begin_facet();
		        builder.add_vertex_to_facet(m_index_counter++);
		        builder.add_vertex_to_facet(m_index_counter++);
		        builder.add_vertex_to_facet(m_index_counter++);
		        builder.add_vertex_to_facet(m_index_counter++);
		        builder.end_facet();
			}
        };

    } // Level_of_detail

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_LOD_1_BUILDER_H