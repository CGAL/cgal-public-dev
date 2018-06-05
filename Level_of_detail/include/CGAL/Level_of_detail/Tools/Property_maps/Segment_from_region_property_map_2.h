#ifndef CGAL_LEVEL_OF_DETAIL_SEGMENT_FROM_REGION_PROPERTY_MAP_2_H
#define CGAL_LEVEL_OF_DETAIL_SEGMENT_FROM_REGION_PROPERTY_MAP_2_H

// CGAL includes.
#include <CGAL/property_map.h>

namespace CGAL {

	namespace Level_of_detail {

		template<typename KeyType, class InputKernel, class Elements, class PointMap>
		class Segment_from_region_property_map_2 {

		public:
			using Kernel 	= InputKernel;
			using Range     = Elements;
			using Point_map = PointMap;
			using ValueType = typename Kernel::Segment_2;
			using Ranges    = std::list<Range>;

			Segment_from_region_property_map_2(const Ranges &input, const Point_map &point_map) :
			m_input(input),
			m_point_map(point_map) 
			{ }

			inline const Point_map& point_map() const {
				return m_point_map;
			}

            using key_type   = KeyType;
            using value_type = ValueType;
            using reference  = const value_type&;
			using category   = boost::lvalue_property_map_tag;
            
            using Self = Point_property_map_2<key_type, value_type, Point_map>;

			value_type operator[](key_type &key) const { 
				return get(this, key);
			}

            friend value_type get(const Self &self, const key_type &key) {
				
            }

		private:
			const Ranges 	&m_input;
			const Point_map &m_point_map;
        };

    } // Level_of_detail

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_SEGMENT_FROM_REGION_PROPERTY_MAP_2_H