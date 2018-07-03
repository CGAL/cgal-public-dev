#ifndef CGAL_LEVEL_OF_DETAIL_BUILDING_WITH_SEGMENT_CONSTRAINTS_PROPERTY_MAP_H
#define CGAL_LEVEL_OF_DETAIL_BUILDING_WITH_SEGMENT_CONSTRAINTS_PROPERTY_MAP_H

// STL includes.
#include <map>
#include <utility>

// CGAL includes.
#include <CGAL/Barycentric_coordinates_2/Segment_coordinates_2.h>

// LOD includes.
#include <CGAL/Level_of_detail/Enumerations.h>

namespace CGAL {

	namespace Level_of_detail {

        namespace LOD = CGAL::Level_of_detail;
        namespace BC  = CGAL::Barycentric_coordinates;

		template<class InputKernel, class InputTriangulation>
		class Building_with_segment_constraints_property_map {
			
        public:
            using Kernel        = InputKernel;
            using Triangulation = InputTriangulation;

            typename Kernel::Compute_squared_distance_2 squared_distance_2;

            using FT        = typename Kernel::FT;
            using Line_2    = typename Kernel::Line_2;
            using Point_2   = typename Kernel::Point_2;
            using Segment_2 = typename Kernel::Segment_2;
            
            using Triangulation_face_handle    = typename Triangulation::Face_handle; 
            using Triangulation_faces_iterator = typename Triangulation::Finite_faces_iterator;
            using All_faces_iterator = typename Triangulation::All_faces_iterator;
            using Triangulation_vertex_handle  = typename Triangulation::Vertex_handle;

            using Face_to_building_map = std::map<Triangulation_face_handle, int>;
            using Visibility_label     = LOD::Visibility_label;

            using Segment_coordinates = std::pair<FT, FT>;

            template<class Elements, class Segment_map>
            Building_with_segment_constraints_property_map(
                const Triangulation &triangulation, 
                const Elements &elements,
                const Segment_map &segment_map,
                const FT constraints_threshold) :
			m_triangulation(triangulation),
            m_constraints_threshold(constraints_threshold),
            m_tolerance(-FT(1) / FT(10)),
            m_building_index(-1),
            m_start_new_building(true) { 

                set_default_building_indices();
                find_buildings(elements, segment_map);
            }

            inline const Face_to_building_map &face_to_building_map() const {
                return m_face_to_building_map;
            }

            inline bool is_valid_face_handle(const Triangulation_face_handle &face_handle, const Face_to_building_map &face_to_building_map) const {
                return face_to_building_map.find(face_handle) != face_to_building_map.end();
            }

            using key_type = Triangulation_face_handle;
            using Self     = Building_with_segment_constraints_property_map<Kernel, Triangulation>;

            friend void put(const Self &self, key_type &face_handle) {
                
                CGAL_precondition(self.is_valid_face_handle(face_handle, self.face_to_building_map()));
                face_handle->info().group_number() = self.face_to_building_map().at(face_handle);
            }

        private:
            const Triangulation &m_triangulation;
            
            const FT m_constraints_threshold;
            const FT m_tolerance;

            int  m_building_index;
            bool m_start_new_building;

            Face_to_building_map m_face_to_building_map;

            void set_default_building_indices() {
                for (All_faces_iterator tf_it = m_triangulation.all_faces_begin(); tf_it != m_triangulation.all_faces_end(); ++tf_it) {
                    
                    const Triangulation_face_handle face_handle = static_cast<Triangulation_face_handle>(tf_it);
                    m_face_to_building_map[face_handle] = -1;
                }
            }

            template<class Elements, class Segment_map>
            void find_buildings(const Elements &elements, const Segment_map &segment_map) {

				for (Triangulation_faces_iterator tf_it = m_triangulation.finite_faces_begin(); tf_it != m_triangulation.finite_faces_end(); ++tf_it) {
                    Triangulation_face_handle face_handle = static_cast<Triangulation_face_handle>(tf_it);
                    
                    flood(elements, segment_map, face_handle);
                    m_start_new_building = true;
				}
            }

            template<class Elements, class Segment_map>
            void flood(const Elements &elements, const Segment_map &segment_map, Triangulation_face_handle &face_handle) {
				
				// If this face is not valid due to some criteria, we do not handle this face and skip it.
				if (is_not_valid_face(face_handle)) return;

				// Check if we should increase the buildings counter.
                set_building_index_to_face(face_handle);

                // Continue propagating.
				for (size_t i = 0; i < 3; ++i) {
				    Triangulation_face_handle face_handle_neighbour = face_handle->neighbor(i);

				    if (!is_constrained_edge(face_handle, i, elements, segment_map)) 
			    		flood(elements, segment_map, face_handle_neighbour);
				}
			}

            bool is_not_valid_face(const Triangulation_face_handle &face_handle) const {

				// This face was already handled before.
				const int building_number = m_face_to_building_map.at(face_handle);
				if (building_number >= 0) return true;

				// This is infinite face.
				if (m_triangulation.is_infinite(face_handle)) return true;

				// This face is outside of any building.
				const Visibility_label visibility_label = face_handle->info().visibility_label();
				if (visibility_label == Visibility_label::OUTSIDE) return true;

				// The face is valid.
				return false;
			}

            void set_building_index_to_face(Triangulation_face_handle &face_handle) {
                
                if (m_start_new_building) {
                    m_start_new_building = false; ++m_building_index;
                }
                m_face_to_building_map[face_handle] = m_building_index;
            }

            template<class Elements, class Segment_map>
            bool is_constrained_edge(const Triangulation_face_handle &face_handle, const size_t vertex_index, const Elements &elements, const Segment_map &segment_map) const {

				const Point_2 &p1 = face_handle->vertex((vertex_index + 1) % 3)->point();
				const Point_2 &p2 = face_handle->vertex((vertex_index + 2) % 3)->point();

				return check_constraint(p1, p2, elements, segment_map);
			}

            template<class Elements, class Segment_map>
            bool check_constraint(const Point_2 &p1, const Point_2 &p2, const Elements &elements, const Segment_map &segment_map) const {
				
				CGAL_precondition(elements.size() > 0);
                using Const_elements_iterator = typename Elements::const_iterator;

                for (Const_elements_iterator ce_it = elements.begin(); ce_it != elements.end(); ++ce_it) {
                    const Segment_2 &segment = get(segment_map, *ce_it);

                    if (is_constrained(p1, p2, segment)) return true;
                }
                return false;
			}

            bool is_constrained(const Point_2 &p1, const Point_2 &p2, const Segment_2 &segment) const {

				const Point_2 &source = segment.source();
				const Point_2 &target = segment.target();

				Line_2 line(source, target);

				const Point_2 pr1 = line.projection(p1);
				const Point_2 pr2 = line.projection(p2);

				const FT squared_threshold = m_constraints_threshold * m_constraints_threshold;

				if (squared_distance_2(p1, pr1) > squared_threshold) return false;
				if (squared_distance_2(p2, pr2) > squared_threshold) return false;
				
				Segment_coordinates bc = CGAL::make_pair(BC::compute_segment_coordinates_2(source, target, p1, Kernel()));
				const bool state1 = bc.first > m_tolerance && bc.second > m_tolerance;

				bc = CGAL::make_pair(BC::compute_segment_coordinates_2(source, target, p2, Kernel()));
				const bool state2 = bc.first > m_tolerance && bc.second > m_tolerance;

				bc = CGAL::make_pair(BC::compute_segment_coordinates_2(p1, p2, source, Kernel()));
				const bool state3 = bc.first > m_tolerance && bc.second > m_tolerance;

				bc = CGAL::make_pair(BC::compute_segment_coordinates_2(p1, p2, target, Kernel()));
				const bool state4 = bc.first > m_tolerance && bc.second > m_tolerance;

				if ( (state1 && state2) || (state3 && state4) ) return true;
				return false;
			}
        };

    } // Level_of_detail

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_BUILDING_WITH_SEGMENT_CONSTRAINTS_PROPERTY_MAP_H
