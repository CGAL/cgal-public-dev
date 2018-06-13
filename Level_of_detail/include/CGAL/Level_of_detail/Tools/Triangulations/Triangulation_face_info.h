#ifndef CGAL_LEVEL_OF_DETAIL_TRIANGULATION_FACE_INFO_H
#define CGAL_LEVEL_OF_DETAIL_TRIANGULATION_FACE_INFO_H

// STL includes.
#include <list>

// LOD includes.
#include <CGAL/Level_of_detail/Enumerations.h>

namespace CGAL {

	namespace Level_of_detail {

        namespace LOD = CGAL::Level_of_detail;

		template<class InputElement>
		class Triangulation_face_info {

		public:
			using Element = InputElement;
            
			using Visibility_label = LOD::Visibility_label;
			using Elements 		   = std::list<Element>;

            Triangulation_face_info() : 
            m_visibility_label(Visibility_label::OUTSIDE),
            m_group_number(-1)
            { }


			// Identifiers.
            inline Visibility_label& visibility_label() {
				return m_visibility_label;
			}

			inline const Visibility_label& visibility_label() const {
				return m_visibility_label;
			}

            inline int& group_number() {
				return m_group_number;
			}

			inline const int& group_number() const {
				return m_group_number;
			}


			// Related elements.
			inline void add_element(const Element &element) {
				m_elements.push_back(element);
			}

			inline Elements& elements() {
				return m_elements;
			}

			inline const Elements& elements() const {
				return m_elements;
			}
            
        private:
            Visibility_label m_visibility_label;
			int              m_group_number;

			Elements m_elements;
        };

    } // Level_of_detail

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_TRIANGULATION_FACE_INFO_H