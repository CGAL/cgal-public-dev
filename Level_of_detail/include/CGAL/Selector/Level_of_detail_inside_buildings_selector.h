#ifndef CGAL_LEVEL_OF_DETAIL_INSIDE_BUILDINGS_SELECTOR_H
#define CGAL_LEVEL_OF_DETAIL_INSIDE_BUILDINGS_SELECTOR_H

#if defined(WIN32) || defined(_WIN32) 
#define PSR "\\"
#else 
#define PSR "/"
#endif

// STL includes.
#include <map>
#include <vector>
#include <unordered_set>

// New CGAL includes.
#include <CGAL/Mylog/Mylog.h>
#include <CGAL/Level_of_detail_enum.h>

namespace CGAL {

	namespace LOD {

		template<class KernelTraits, class InputContainer, class InputCDT, class InputBuildings>
		class Level_of_detail_inside_buildings_selector {

        public:
            typedef KernelTraits   Kernel;
            typedef InputContainer Input;
            typedef InputCDT       CDT;
            typedef InputBuildings Buildings;
			
            using Index   = int;
			using Indices = std::vector<Index>;

            using Vertex_handle = typename CDT::Vertex_handle;
            using Face_handle   = typename CDT::Face_handle;

            using FT      = typename Kernel::FT;
            using Point_3 = typename Kernel::Point_3;

            using Point_2    = typename Kernel::Point_2;
            using Triangle_2 = typename Kernel::Triangle_2;

            using Faces = std::vector<Face_handle>;

            using Building          = CGAL::LOD::Building<FT, Vertex_handle, Face_handle, Point_3>;
            using Building_iterator = typename Buildings::iterator;

            using Log = CGAL::LOD::Mylog;

            Level_of_detail_inside_buildings_selector(const Input &input, const CDT &cdt, const Indices &indices) : 
            m_input(input), m_cdt(cdt), m_indices(indices),
            m_silent(false), m_height_threshold(FT(1) / FT(1000000)) { }

            void add_indices(Buildings &buildings) const {
                if (buildings.size() == 0) return;

                for (Building_iterator bit = buildings.begin(); bit != buildings.end(); ++bit) {
                    Building &building = (*bit).second;

                    if (is_valid_building(buildings, building)) {
                        add_indices_to_building(building);
                        
                        if (building.interior_indices.size() < 3) 
                            building.is_valid = false;
                    }
                }

                if (!m_silent) {
                    Log exporter; exporter.export_points_inside_buildings(buildings, m_input, "tmp" + std::string(PSR) + "inside_buildings_points");
                }
            }

            void make_silent(const bool new_state) {
                m_silent = new_state;
            }

        private:
            const Input   &m_input;
            const CDT     &m_cdt;
            const Indices &m_indices;

            bool  m_silent;
            const FT m_height_threshold;

			bool is_valid_building(const Buildings &buildings, Building &building) const {

				const FT height = building.height;
				if (height < m_height_threshold) return false;

				const auto &faces = building.faces;
				if (faces.size() < 2) {
				
					for (std::unordered_set<int>::const_iterator nit = building.neighbours.begin(); nit != building.neighbours.end(); ++nit) {
						if (is_valid_local_building(buildings.at(*nit))) {

							building.height = buildings.at(*nit).height;
							building.color  = buildings.at(*nit).color;

							return true;
						}
					}
					return false;
				}
				return true;
			}

			bool is_valid_local_building(const Building &building) const {

				const FT height = building.height;
				if (height < m_height_threshold) return false;

				const auto &faces = building.faces;
				if (faces.size() < 2) return false;

				return true;
			}

            void add_indices_to_building(Building &building) const {
                const Faces &building_faces = building.faces;
                
                for (size_t i = 0; i < building_faces.size(); ++i)
                    add_face_indices(building_faces[i], building);
            }

            void add_face_indices(const Face_handle &fh, Building &building) const {

                for (size_t i = 0; i < m_indices.size(); ++i) {
                    const Index index = m_indices[i];

                    const Point_3 &query = m_input.point(index);
                    if (belongs_to_face(query, fh)) building.interior_indices.push_back(index); 
                }
            }

            bool belongs_to_face(const Point_3 &query, const Face_handle &fh) const {

				const Point_2 &p1 = fh->vertex(0)->point();
				const Point_2 &p2 = fh->vertex(1)->point();
				const Point_2 &p3 = fh->vertex(2)->point();

				const Triangle_2 triangle = Triangle_2(p1, p2, p3);
				const Point_2 new_query   = Point_2(query.x(), query.y());

				if (triangle.has_on_bounded_side(new_query) || triangle.has_on_boundary(new_query)) return true;
				return false;
			}
        };
    }
}

#endif // CGAL_LEVEL_OF_DETAIL_INSIDE_BUILDINGS_SELECTOR_H