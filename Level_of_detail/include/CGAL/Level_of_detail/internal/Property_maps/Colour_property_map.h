#ifndef CGAL_LEVEL_OF_DETAIL_COLOUR_PROPERTY_MAP_H
#define CGAL_LEVEL_OF_DETAIL_COLOUR_PROPERTY_MAP_H

// STL includes.
#include <time.h>
#include <stdio.h>
#include <stdlib.h>

// CGAL includes.
#include <CGAL/IO/Color.h>
#include <CGAL/property_map.h>

// LOD includes.
#include <CGAL/Level_of_detail/Enumerations.h>

namespace CGAL {

	namespace Level_of_detail {

		namespace LOD = CGAL::Level_of_detail;

		class Colour_property_map {

		public:
            using Colour    = CGAL::Color;
            using ValueType = Colour;

			using Colour_map_type = LOD::Colour_map_type;

			Colour_property_map(const Colour_map_type colour_map_type) :
			m_colour_map_type(colour_map_type) { 
				
				srand(time(NULL));
			}

            using value_type = ValueType;
            using reference  = const value_type&;
			using category   = boost::lvalue_property_map_tag;
            
            using Self = Colour_property_map;

			template<typename key_type>
			value_type operator[](key_type &) const { }

			Colour_map_type colour_map_type() const {
				return m_colour_map_type;
			}

			template<typename key_type>
            friend value_type get(const Self &self, const key_type &) {
				
				switch (self.colour_map_type()) {
					
					case Colour_map_type::RANDOM:
						return self.generate_random_colour();

					case Colour_map_type::WHITE:
						return self.generate_white_colour();

					case Colour_map_type::BLACK:
						return self.generate_black_colour();

					case Colour_map_type::GROUND_DEFAULT:
						return self.generate_ground_default_colour();

					case Colour_map_type::WALL_DEFAULT:
						return self.generate_wall_default_colour();

					case Colour_map_type::ROOF_DEFAULT:
						return self.generate_roof_default_colour();
					
					default:
						return self.generate_random_colour();
				}
            }

            Colour generate_random_colour() const {

				const int r = rand() % 255;
				const int g = rand() % 255;
				const int b = rand() % 255;

				return Colour(r, g, b);
			}

			Colour generate_white_colour() const {
				return Colour(255, 255, 255);
			}

			Colour generate_black_colour() const {
				return Colour(0, 0, 0);
			}

			Colour generate_ground_default_colour() const {
				return Colour(204, 191, 144);
			}

			Colour generate_wall_default_colour() const {
				return Colour(114, 117, 132);
			}

			Colour generate_roof_default_colour() const {
				return Colour(145, 84, 93);
			}

        private:
			const Colour_map_type m_colour_map_type;
		};

    } // Level_of_detail

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_COLOUR_PROPERTY_MAP_H