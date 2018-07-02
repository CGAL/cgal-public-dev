#ifndef CGAL_LEVEL_OF_DETAIL_DEREFERENCE_PROPERTY_MAP_H
#define CGAL_LEVEL_OF_DETAIL_DEREFERENCE_PROPERTY_MAP_H

// CGAL includes.
#include <CGAL/property_map.h>

namespace CGAL {

	namespace Level_of_detail {

		template<typename KeyType, class PropertyMap>
		class Dereference_property_map {

		public:
			using Property_map = PropertyMap;

			Dereference_property_map(const Property_map& property_map) :
        m_property_map(property_map) 
			{ }

			inline const Property_map& property_map() const {
				return m_property_map;
			}

      using key_type   = KeyType;
      using value_type = typename Property_map::value_type;
      using reference  = const value_type&;
			using category   = boost::lvalue_property_map_tag;
            
      using Self = Dereference_property_map<key_type, Property_map>;

			reference operator[](key_type &key) const { 
				return get(this, key);
			}

      friend reference get(const Self &self, const key_type &key) {
				
        return get(self.property_map(), *key);
      }

    private:
      Property_map m_property_map;
		};

  } // Level_of_detail

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_DEREFERENCE_PROPERTY_MAP_H
