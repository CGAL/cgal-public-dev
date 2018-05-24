#ifndef CGAL_LEVEL_OF_DETAIL_DATA_STRUCTURE_H
#define CGAL_LEVEL_OF_DETAIL_DATA_STRUCTURE_H

// STL indludes.
#include <list>

namespace CGAL {

	namespace Level_of_detail {

		template<class InputKernel, class InputData, class InputPointMap, class InputNormalMap, class InputLabelMap>
		struct Level_of_detail_data_structure {

		public:			
            using Kernel     = InputKernel;
			using Input      = InputData;
            using Point_map  = InputPointMap;
			using Normal_map = InputNormalMap;
			using Label_map  = InputLabelMap;

            using Point_index   = typename Input::const_iterator;
			using Point_indices = std::list<Point_index>;

            Level_of_detail_data_structure(const Input &input, const Point_map &point_map, const Normal_map &normal_map, const Label_map &label_map) :
            m_input(input),
            m_point_map(point_map),
            m_normal_map(normal_map),
            m_label_map(label_map)
            { }

            // Input.
            inline const Input& input() const {
                return m_input;
            }

            inline const Point_map& point_map() const {
                return m_point_map;
            }

            inline const Normal_map& normal_map() const {
                return m_normal_map;
            }

            inline const Label_map& label_map() const {
                return m_label_map;
            }

            // Indices.
            inline Point_indices& ground_indices() {
                return m_ground_indices;
            }

            inline const Point_indices& ground_indices() const {
                return m_ground_indices;
            }

            inline Point_indices& building_boundary_indices() {
                return m_building_boundary_indices;
            }

            inline const Point_indices& building_boundary_indices() const {
                return m_building_boundary_indices;
            }

            inline Point_indices& building_interior_indices() {
                return m_building_interior_indices;
            }

            inline const Point_indices& building_interior_indices() const {
                return m_building_interior_indices;
            }

        private:
            const Input      &m_input;
            const Point_map  &m_point_map;
            const Normal_map &m_normal_map;
            const Label_map  &m_label_map;

            Point_indices m_ground_indices;
            Point_indices m_building_boundary_indices;
            Point_indices m_building_interior_indices;
        };
    
    } // Level_of_detail

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_DATA_STRUCTURE_H