#ifndef CGAL_LEVEL_OF_DETAIL_LOD2_H
#define CGAL_LEVEL_OF_DETAIL_LOD2_H

// STL includes.
#include <map>
#include <string>
#include <vector>
#include <cassert>
#include <iostream>

// CGAL includes.
#include <CGAL/utils.h>
#include <CGAL/Random.h>
#include <CGAL/IO/Color.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/number_utils.h>
#include <CGAL/utils_classes.h>
#include <CGAL/Unique_hash_map.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>

// New CGAL includes.
#include <CGAL/Mylog/Mylog.h>
#include <CGAL/Level_of_detail_enum.h>
#include <CGAL/Utils/Level_of_detail_utils.h>

namespace CGAL {

	namespace LOD {

		template<class KernelTraits, class InputHDS, class InputBuilding, class InputBuildings, class FacetHandle>
		class LOD2_builder : public Modifier_base<InputHDS> {
		
		public:
			typedef KernelTraits   Kernel;
    		typedef InputHDS 	   HDS;
    		typedef InputBuilding  Building;
			typedef InputBuildings Buildings;
    		typedef FacetHandle    Color_facet_handle;
    		typedef FacetHandle    Facet_handle;

    		using FT 	  = typename Kernel::FT;
			using Point_2 = typename Kernel::Point_2;
			using Point_3 = typename Kernel::Point_3;

			using Boundary = typename Building::Boundary;			
			using Roof 	   = typename Building::Roof;
			using Roofs    = typename Building::Roofs;
			
			using Roof_boundary = typename Roof::Roof_boundary;
			using Polygon 		= Roof_boundary;

			using Color 	   = CGAL::Color;
			using Facet_colors = std::map<Color_facet_handle, Color>;
			using Builder 	   = CGAL::Polyhedron_incremental_builder_3<HDS>;
    		
			using Building_iterator = typename Buildings::const_iterator;
    		using Ground 			= std::vector<Point_3>;

			using Log 		   = CGAL::LOD::Mylog;
			using Simple_utils = Level_of_detail_utils_simple<Kernel>;

			LOD2_builder(const Buildings &buildings, const Ground &ground, const FT ground_height, Facet_colors &facet_colors) : 
    		m_buildings(buildings), 
			m_ground(ground),
			m_ground_height(ground_height),
			m_facet_colors(facet_colors),
    		m_index_counter(0) { }

			void operator()(HDS &hds) {
				
				Builder builder(hds, false);

				const size_t expected_num_vertices = estimate_number_of_vertices();
				const size_t expected_num_facets   = estimate_number_of_facets();

				m_index_counter = 0;
				builder.begin_surface(expected_num_vertices, expected_num_facets);

				for (Building_iterator bit = m_buildings.begin(); bit != m_buildings.end(); ++bit) {
				
					const Building &building = (*bit).second;
					if (!building.is_valid) continue;

					add_new_building(building, builder);
				}

				add_ground(builder);
				builder.end_surface();
			}

		private:
			const Buildings &m_buildings;
			const Ground    &m_ground;
			
			const FT m_ground_height;
			Facet_colors &m_facet_colors;
			
			size_t 		 m_index_counter;
			Simple_utils m_simple_utils;

			inline size_t estimate_number_of_vertices() {
				return m_buildings.size() * 8 + 4;
			}

			inline size_t estimate_number_of_facets() {
				return m_buildings.size() * 6 + 1;
			}

			void add_new_building(const Building &building, Builder &builder) {

				const Roofs &roofs = building.roofs;
				const Color &color = building.color;
				
				add_structure_from_roofs(roofs, color, builder);
			}

			void add_structure_from_roofs(const Roofs &roofs, const Color &color, Builder &builder) {
				
				const size_t num_roofs = roofs.size();
				assert(num_roofs > 0);

				for (size_t i = 0; i < num_roofs; ++i) {
					
					const Roof &roof = roofs[i];
					add_structure_from_roof(roof, color, builder);
				}
			}

			void add_structure_from_roof(const Roof &roof, const Color &color, Builder &builder) {

				const Roof_boundary &boundary = roof.boundary;
				assert(boundary.size() > 2);

				 add_roof(boundary, color, builder);
				add_walls(boundary, color, builder);
			}

			void add_roof(const Polygon &vertices, const Color &color, Builder &builder) {

				for (size_t i = 0; i < vertices.size(); ++i)
					builder.add_vertex(vertices[i]);

				const Color_facet_handle cfh = builder.begin_facet();
				for (size_t i = 0; i < vertices.size(); ++i)
					builder.add_vertex_to_facet(m_index_counter++);
				builder.end_facet();
				
		        m_facet_colors[cfh] = color;
			}

			void add_walls(const Polygon &vertices, const Color &color, Builder &builder) {
				const size_t n = vertices.size();

				for (size_t i = 0; i < n; ++i) {
					const size_t ip = (i + 1) % n;

					const Point_3 &a = vertices[i];
					const Point_3 &b = vertices[ip];

					const Point_3 c = Point_3(b.x(), b.y(), m_ground_height);
					const Point_3 d = Point_3(a.x(), a.y(), m_ground_height);
					
					add_quad(a, b, c, d, color, builder);
				}		
			}

			void add_ground(Builder &builder) {
				assert(!m_ground.empty());

				const size_t num_vertices = m_ground.size();
				assert(num_vertices == 4);

				const Point_3 &a = m_ground[3];
				const Point_3 &b = m_ground[2];
				const Point_3 &c = m_ground[1];
				const Point_3 &d = m_ground[0];

				const Point_3 p1 = Point_3(a.x(), a.y(), a.z() + m_ground_height);
				const Point_3 p2 = Point_3(b.x(), b.y(), b.z() + m_ground_height);
				const Point_3 p3 = Point_3(c.x(), c.y(), c.z() + m_ground_height);
				const Point_3 p4 = Point_3(d.x(), d.y(), d.z() + m_ground_height);

				const Color color(169, 169, 169);
				add_quad(p1, p2, p3, p4, color, builder);
			}

			void add_quad(const Point_3 &a, const Point_3 &b, const Point_3 &c, const Point_3 &d, const Color &color, Builder &builder) {

		        builder.add_vertex(a);
		        builder.add_vertex(b);
		        builder.add_vertex(c);
		        builder.add_vertex(d);

		        const Color_facet_handle cfh = builder.begin_facet();

		        builder.add_vertex_to_facet(m_index_counter++);
		        builder.add_vertex_to_facet(m_index_counter++);
		        builder.add_vertex_to_facet(m_index_counter++);
		        builder.add_vertex_to_facet(m_index_counter++);

		        builder.end_facet();
		        m_facet_colors[cfh] = color;
			}
        };

        template<class KernelTraits, class InputBuilding, class InputBuildings, class OutputMesh>
		class Level_of_detail_lod2 {

		public:
			typedef KernelTraits   Kernel;
			typedef InputBuilding  Building;
			typedef InputBuildings Buildings;
			typedef OutputMesh 	   Mesh;

			using FT  = typename Kernel::FT;
            using HDS = typename Mesh::HalfedgeDS;

            using Mesh_facet_handle = typename Mesh::Facet_const_handle;
            using Mesh_builder 		= LOD2_builder<Kernel, HDS, Building, Buildings, Mesh_facet_handle>;
			using Mesh_facet_colors = typename Mesh_builder::Facet_colors;
			using Ground 			= typename Mesh_builder::Ground;

            Level_of_detail_lod2(const Buildings &buildings, const Ground &ground, const FT ground_height, Mesh_facet_colors &mesh_facet_colors) :
            m_builder(buildings, ground, ground_height, mesh_facet_colors) { }

            inline void reconstruct(Mesh &mesh) {
				mesh.delegate(m_builder);
            }

        private:
			Mesh_builder m_builder;
        };
    }
}

#endif // CGAL_LEVEL_OF_DETAIL_LOD2_H