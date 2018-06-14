#ifndef CGAL_LEVEL_OF_DETAIL_VISIBILITY_FROM_FACETS_COLOUR_PROPERTY_MAP_H
#define CGAL_LEVEL_OF_DETAIL_VISIBILITY_FROM_FACETS_COLOUR_PROPERTY_MAP_H

// CGAL includes.
#include <CGAL/IO/Color.h>
#include <CGAL/property_map.h>

// LOD includes.
#include <CGAL/Level_of_detail/Enumerations.h>

namespace CGAL {

	namespace Level_of_detail {

		namespace LOD = CGAL::Level_of_detail;

		class Visibility_from_facets_colour_property_map {

		public:
            using Colour    = CGAL::Color;
            using ValueType = Colour;

            using value_type = ValueType;
            using reference  = const value_type&;
			using category   = boost::lvalue_property_map_tag;
            
            using Self = Visibility_from_facets_colour_property_map;
			using Visibility_label = LOD::Visibility_label;

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
				
				const Visibility_label visibility_label = facet.info().visibility_label();
				switch (visibility_label) {

					case Visibility_label::OUTSIDE:
						return Colour(255, 55, 55); // red

					case Visibility_label::INSIDE:
						return Colour(55, 255, 55); // green

					default:
						return Colour(255, 215, 0); // yellow
				}
			}
		};

    } // Level_of_detail

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_VISIBILITY_FROM_FACETS_COLOUR_PROPERTY_MAP_H