#ifndef CGAL_LEVEL_OF_DETAIL_ESTIMATED_NORMAL_PROPERTY_MAP_2_H
#define CGAL_LEVEL_OF_DETAIL_ESTIMATED_NORMAL_PROPERTY_MAP_2_H

// CGAL includes.
#include <CGAL/property_map.h>
#include <CGAL/number_utils.h>

namespace CGAL {

	namespace Level_of_detail {

		template<class InputKernel, class LinesTraits_2>
		class Estimated_normal_property_map_2 {

		public:
            using Kernel         = InputKernel;
            using Lines_traits_2 = LinesTraits_2;

            using FT        = typename Kernel::FT;
            using Vector_2  = typename Kernel::Vector_2;
            using Line_2    = typename Kernel::Line_2;
            using ValueType = Vector_2;

            using Point_identifier = typename Lines_traits_2::Point_identifier;
            
            using KeyType = Point_identifier;
            using Normals = std::map<KeyType, Vector_2>;

            using Elements  = typename Lines_traits_2::Elements;
            using Point_map = typename Lines_traits_2::Point_map;
            using Lines_2   = typename Lines_traits_2::Lines_2;

            using Elements_iterator = typename Elements::const_iterator;

			Estimated_normal_property_map_2(const Elements &elements, const Point_map &point_map, const Lines_2 &lines_2) :
			m_elements(elements),
            m_point_map(point_map),
            m_lines_2(lines_2) { 

                estimate_normals();
            }

            const Normals &normals() const {
                return m_normals;
            }

            using key_type   = KeyType;
            using value_type = ValueType;
            using reference  = const value_type&;
            using category   = boost::lvalue_property_map_tag;
            
            using Self = Estimated_normal_property_map_2<Kernel, Lines_traits_2>;

			reference operator[](key_type &key) const { 
				return get(this, key);
			}

            friend reference get(const Self &self, const key_type &key) {
				return self.normals().at(key);
            }

        private:
            const Elements  &m_elements;
            const Point_map &m_point_map;
            const Lines_2   &m_lines_2;

            Normals m_normals;

            void estimate_normals() {
                
                m_normals.clear();
                CGAL_precondition(m_elements.size() > 0);
                
                for (Elements_iterator element = m_elements.begin(); element != m_elements.end(); ++element) 
                    estimate_normal(*element);
            }

            void estimate_normal(const Point_identifier &point_id) {

                const Line_2 &line    = m_lines_2.at(point_id);
				const Vector_2 vector = line.to_vector();
                
				Vector_2 normal       = vector.perpendicular(CGAL::COUNTERCLOCKWISE);
				const FT length       = static_cast<FT>(CGAL::sqrt(CGAL::to_double(normal * normal)));

				CGAL_precondition(length != FT(0));
				normal /= length;
                m_normals[point_id] = normal;
            }
		};

    } // Level_of_detail

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_ESTIMATED_NORMAL_PROPERTY_MAP_2_H