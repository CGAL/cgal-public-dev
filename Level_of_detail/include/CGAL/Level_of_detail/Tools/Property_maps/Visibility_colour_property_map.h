#ifndef CGAL_LEVEL_OF_DETAIL_VISIBILITY_COLOUR_PROPERTY_MAP_H
#define CGAL_LEVEL_OF_DETAIL_VISIBILITY_COLOUR_PROPERTY_MAP_H

// CGAL includes.
#include <CGAL/IO/Color.h>
#include <CGAL/property_map.h>

namespace CGAL {

	namespace Level_of_detail {

		class Visibility_colour_property_map {

		public:
            using Colour    = CGAL::Color;
            using ValueType = Colour;

            using value_type = ValueType;
            using reference  = const value_type&;
			using category   = boost::lvalue_property_map_tag;
            
            using Self = Visibility_colour_property_map;

			template<typename key_type>
			value_type operator[](key_type &key) const { 
				return get(this, key);
			}

			template<typename key_type>
            friend value_type get(const Self &self, const key_type &key) {
				return self.generate_visibility_colour(key);
            }

			template<class Facet>
			Colour generate_visibility_colour(const Facet &facet) const {
				
				const size_t value = static_cast<size_t>(CGAL::to_double(facet.building_interior()));
				switch (value) {

					case 0:
						return Colour(255, 55, 55); // red

					case 1:
						return Colour(55, 255, 55); // green

					default:
						return Colour(255, 215, 0); // yellow
				}
			}
		};

    } // Level_of_detail

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_VISIBILITY_COLOUR_PROPERTY_MAP_H