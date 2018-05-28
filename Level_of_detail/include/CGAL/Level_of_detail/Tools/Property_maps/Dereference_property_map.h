#ifndef CGAL_LEVEL_OF_DETAIL_DEREFERENCE_PROPERTY_MAP_H
#define CGAL_LEVEL_OF_DETAIL_DEREFERENCE_PROPERTY_MAP_H

namespace CGAL {

	namespace Level_of_detail {

		template<typename KeyType, class PropertyMap>
		class Dereference_property_map {

		public:
			using Property_map = PropertyMap;

			Dereference_property_map(const Property_map &property_map) :
			m_property_map(property_map) 
			{ }

			inline const Property_map& property_map() const {
				return m_property_map;
			}

		private:
            using key_type   = KeyType;
            using value_type = typename Property_map::value_type;
            using reference  = const value_type&;
            
            using Self = Dereference_property_map<key_type, Property_map>;

            friend reference get(const Self &self, const key_type &key) {
				
				const Property_map &property_map = self.property_map();
                return get(property_map, *key);
            }

            friend void put(const Self &, key_type &, const value_type &) { }
			const Property_map &m_property_map;
		};

    } // Level_of_detail

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_DEREFERENCE_PROPERTY_MAP_H