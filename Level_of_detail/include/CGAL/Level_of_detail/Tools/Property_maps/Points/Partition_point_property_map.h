#ifndef CGAL_LEVEL_OF_DETAIL_PARTITION_POINT_PROPERTY_MAP_H
#define CGAL_LEVEL_OF_DETAIL_PARTITION_POINT_PROPERTY_MAP_H

// CGAL includes.
#include <CGAL/property_map.h>

namespace CGAL {

	namespace Level_of_detail {

		template<typename KeyType, typename ValueType>
		class Partition_point_property_map {

		public:
            using key_type   = KeyType;
            using value_type = ValueType;
            using reference  = const value_type&;
			using category   = boost::lvalue_property_map_tag;
            
            using Self = Partition_point_property_map<key_type, value_type>;

			value_type operator[](key_type &key) const { 
				return get(this, key);
			}

            friend value_type get(const Self &self, const key_type &key) {
				return value_type(key.x(), key.y(), 0);
            }
		};

    } // Level_of_detail

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_PARTITION_POINT_PROPERTY_MAP_H