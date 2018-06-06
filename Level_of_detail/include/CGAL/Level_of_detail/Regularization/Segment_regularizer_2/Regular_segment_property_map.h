#ifndef CGAL_LEVEL_OF_DETAIL_REGULAR_SEGMENT_PROPERTY_MAP_H
#define CGAL_LEVEL_OF_DETAIL_REGULAR_SEGMENT_PROPERTY_MAP_H

namespace CGAL {

	namespace Level_of_detail {

		template<typename KeyType, typename ValueType>
		class Regular_segment_property_map {
			
		private:
            typedef KeyType          key_type;
            typedef ValueType        value_type;
            typedef const ValueType& reference;

            typedef Regular_segment_property_map<KeyType, ValueType> Self;

            friend reference get(const Self&, const key_type* key) {
                return key->get();
            }

            friend void put(const Self&, key_type &key, const value_type &value) {
                key.get() = value;
            }
		};

	} // Level_of_detail

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_REGULAR_SEGMENT_PROPERTY_MAP_H