#ifndef CGAL_LEVEL_OF_DETAIL_DATA_STRUCTURE_H
#define CGAL_LEVEL_OF_DETAIL_DATA_STRUCTURE_H

// STL indludes.
#include <list>

namespace CGAL {

	namespace Level_of_detail {

		template<class InputKernel, class InputRange, class InputPointMap>
		struct Data_structure {

		public:			
            using Kernel      = InputKernel;
			using Input_range = InputRange;
            using Point_map   = InputPointMap;

            using Point_index   = typename Input_range::const_iterator;
			using Point_indices = std::list<Point_index>;

            Data_structure(const Input_range &input_range, const Point_map &point_map) :
            m_input_range(input_range),
            m_point_map(point_map)
            { }

            // Input.
            inline const Input_range& input_range() const {
                return m_input_range;
            }

            inline const Point_map& point_map() const {
                return m_point_map;
            }

            // Indices.
            inline Point_indices& ground_points() {
                return m_ground_indices;
            }

            inline const Point_indices& ground_points() const {
                return m_ground_indices;
            }

            inline Point_indices& building_boundary_points() {
                return m_building_boundary_indices;
            }

            inline const Point_indices& building_boundary_points() const {
                return m_building_boundary_indices;
            }

            inline Point_indices& building_interior_points() {
                return m_building_interior_indices;
            }

            inline const Point_indices& building_interior_points() const {
                return m_building_interior_indices;
            }

        private:
            const Input_range &m_input_range;
            const Point_map   &m_point_map;

            Point_indices m_ground_indices;
            Point_indices m_building_boundary_indices;
            Point_indices m_building_interior_indices;
        };
    
    } // Level_of_detail

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_DATA_STRUCTURE_H