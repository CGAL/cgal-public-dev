#ifndef CGAL_LEVEL_OF_DETAIL_POINT_PROPERTY_MAP_2_H
#define CGAL_LEVEL_OF_DETAIL_POINT_PROPERTY_MAP_2_H

namespace CGAL {

	namespace Level_of_detail {

		template<typename KeyType, typename ValueType, class PointMap>
		class Point_property_map_2 {

		public:
			using Point_map = PointMap;

			Point_property_map_2(const Point_map &point_map) :
			m_point_map(point_map) 
			{ }

			inline const Point_map& point_map() const {
				return m_point_map;
			}

		private:
            using key_type   = KeyType;
            using value_type = ValueType;
            using reference  = const value_type&;
            
            using Self = Point_property_map_2<key_type, value_type, Point_map>;

            friend value_type get(const Self &self, const key_type &key) {
				
				const Point_map &point_map = self.point_map();
                const auto& point = get(point_map, *key);

				return value_type(point.x(), point.y());
            }

            friend void put(const Self &, key_type &, const value_type &) { }
			const Point_map &m_point_map;
		};

    } // Level_of_detail

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_POINT_PROPERTY_MAP_2_H