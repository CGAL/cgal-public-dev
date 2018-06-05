#ifndef CGAL_LEVEL_OF_DETAIL_DATA_STRUCTURE_H
#define CGAL_LEVEL_OF_DETAIL_DATA_STRUCTURE_H

// STL indludes.
#include <list>
#include <vector>

namespace CGAL {

	namespace Level_of_detail {

		template<class InputKernel, class InputRange, class InputPointMap>
		struct Data_structure {

		public:			
            using Kernel      = InputKernel;
			using Input_range = InputRange;
            using Point_map   = InputPointMap;

            using Point_identifier  = typename Input_range::const_iterator;
			using Point_identifiers = std::vector<Point_identifier>;

            using Point_3   = typename Kernel::Point_3;
            using Plane_3   = typename Kernel::Plane_3;
            using Segment_2 = typename Kernel::Segment_2;

            using Polygon_3 = std::vector<Point_3>;
            
            using Detected_regions     = std::list<Point_identifiers>;
            using Regularized_segments = std::list<Segment_2>;

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

            // Points.
            inline Point_identifiers& ground_points() {
                return m_ground_point_identifiers;
            }

            inline const Point_identifiers& ground_points() const {
                return m_ground_point_identifiers;
            }

            inline Point_identifiers& building_boundary_points() {
                return m_building_boundary_point_identifiers;
            }

            inline const Point_identifiers& building_boundary_points() const {
                return m_building_boundary_point_identifiers;
            }

            inline Point_identifiers& building_interior_points() {
                return m_building_interior_point_identifiers;
            }

            inline const Point_identifiers& building_interior_points() const {
                return m_building_interior_point_identifiers;
            }

            inline Point_identifiers& filtered_building_boundary_points() {
                return m_filtered_building_boundary_point_identifiers;
            }

            inline const Point_identifiers& filtered_building_boundary_points() const {
                return m_filtered_building_boundary_point_identifiers;
            }

            inline Point_identifiers& simplified_building_boundary_points() {
                return m_simplified_building_boundary_point_identifiers;
            }

            inline const Point_identifiers& simplified_building_boundary_points() const {
                return m_simplified_building_boundary_point_identifiers;
            }

            // Ground.
            inline Plane_3& ground_plane() {
                return m_ground_plane;
            }

            inline const Plane_3& ground_plane() const {
                return m_ground_plane;
            }

            inline Polygon_3& ground_bounding_box() {
                return m_ground_bounding_box;
            }

            inline const Polygon_3& ground_bounding_box() const {
                return m_ground_bounding_box;
            }

            // Regions.
            inline Detected_regions& detected_2d_regions() {
                return m_detected_2d_regions;
            }

            inline const Detected_regions& detected_2d_regions() const {
                return m_detected_2d_regions;
            }

            // Regularized segments.
            inline Regularized_segments& regularized_segments() {
                return m_regularized_segments;
            }

            inline const Regularized_segments& regularized_segments() const {
                return m_regularized_segments;
            }

        private:
            const Input_range &m_input_range;
            const Point_map   &m_point_map;

            Point_identifiers m_ground_point_identifiers;
            Point_identifiers m_building_boundary_point_identifiers;
            Point_identifiers m_building_interior_point_identifiers;

            Plane_3   m_ground_plane;
            Polygon_3 m_ground_bounding_box;

            Point_identifiers m_filtered_building_boundary_point_identifiers;
            Point_identifiers m_simplified_building_boundary_point_identifiers;

            Detected_regions     m_detected_2d_regions;
            Regularized_segments m_regularized_segments;
        };
    
    } // Level_of_detail

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_DATA_STRUCTURE_H