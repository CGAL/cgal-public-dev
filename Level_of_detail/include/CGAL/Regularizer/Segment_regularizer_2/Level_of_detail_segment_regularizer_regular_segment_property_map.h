#ifndef CGAL_LEVEL_OF_DETAIL_SEGMENT_REGULARIZER_REGULAR_SEGMENT_PROPERTY_MAP_H
#define CGAL_LEVEL_OF_DETAIL_SEGMENT_REGULARIZER_REGULAR_SEGMENT_PROPERTY_MAP_H

namespace CGAL {

	namespace LOD {

		template<typename KeyType, typename ValueType>
		class Level_of_detail_segment_regularizer_regular_segment_property_map {
			
		private:
            typedef KeyType          key_type;
            typedef ValueType        value_type;
            typedef const ValueType& reference;

            typedef Level_of_detail_segment_regularizer_regular_segment_property_map<KeyType, ValueType> Self;

            friend reference get(const Self&, const key_type& key) {
                return key.get();
            }

            friend void put(const Self&, key_type &key, const value_type &value) {
                key.get() = value;
            }
		};

	} // namespace LOD

} // namespace CGAL

#endif // CGAL_LEVEL_OF_DETAIL_SEGMENT_REGULARIZER_REGULAR_SEGMENT_PROPERTY_MAP_H