#ifndef CGAL_LEVEL_OF_DETAIL_BUILDINGS_OUTLINER_H
#define CGAL_LEVEL_OF_DETAIL_BUILDINGS_OUTLINER_H

// LOD includes.
#include <CGAL/Level_of_detail/Enumerations.h>

namespace CGAL {

	namespace Level_of_detail {

        namespace LOD = CGAL::Level_of_detail;

        template<class InputKernel, class InputBuilding>
        class Buildings_outliner {
        
        public:
            using Kernel   = InputKernel;
            using Building = InputBuilding;

            using Point_2 = typename Kernel::Point_2;

            using Face_handle           = typename Building::Floor_face_handle;
            using Face_handles_iterator = typename Building::Const_floor_face_handles_iterator;

            using Visibility_label = LOD::Visibility_label;

            template<class Triangulation, class Buildings>
            void find_walls(const Triangulation &triangulation, Buildings &buildings) const {
                using Buildings_iterator = typename Buildings::iterator;

                CGAL_precondition(buildings.size() > 0);
                for (Buildings_iterator bu_it = buildings.begin(); bu_it != buildings.end(); ++bu_it) {
                    
                    Building &building = *bu_it;
					find_building_boundaries(triangulation, building);
				}
            }

            template<class Triangulation>
            void find_building_boundaries(const Triangulation &triangulation, Building &building) const {
                for (Face_handles_iterator fh_it = building.floor_face_handles().begin(); fh_it != building.floor_face_handles().end(); ++fh_it) {
                    
                    const Face_handle &face_handle = *fh_it;
                    add_floor_face_edges(triangulation, face_handle, building);
                }
            }

            template<class Triangulation>
            void add_floor_face_edges(const Triangulation &triangulation, const Face_handle &face_handle, Building &building) const {

				for (size_t i = 0; i < 3; ++i) {
					const Face_handle &face_handle_neighbour = face_handle->neighbor(i);

					if (is_boundary_edge_of_building(triangulation, face_handle, face_handle_neighbour))
						add_floor_face_edge(i, face_handle, building);
				}
            }

            template<class Triangulation>
			bool is_boundary_edge_of_building(const Triangulation &triangulation, const Face_handle &face_handle, const Face_handle &face_handle_neighbour) const {

				// This is infinite face.
				if (triangulation.is_infinite(face_handle_neighbour)) return true;

                // This not a building at all.
                const Visibility_label visibility_label = face_handle_neighbour->info().visibility_label();
				if (visibility_label == Visibility_label::OUTSIDE) return true;

                // This is another building.
				const int building_index           = face_handle->info().group_number();
				const int building_index_neighbour = face_handle_neighbour->info().group_number();

				if (building_index_neighbour >= 0 && building_index_neighbour != building_index) 
                    return true;

                // This not a boundary edge.
				return false;
			}

            void add_floor_face_edge(const size_t vertex_index, const Face_handle &face_handle, Building &building) const {

				const Point_2 &p1 = face_handle->vertex((vertex_index + 1) % 3)->point();
				const Point_2 &p2 = face_handle->vertex((vertex_index + 2) % 3)->point();

				building.add_floor_edge(p1, p2);
            }
        };

    } // Level_of_detail

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_BUILDINGS_OUTLINER_H