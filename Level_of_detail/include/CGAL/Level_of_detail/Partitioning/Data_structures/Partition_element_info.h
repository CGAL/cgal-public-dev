#ifndef CGAL_LEVEL_OF_DETAIL_PARTITION_ELEMENT_INFO_H
#define CGAL_LEVEL_OF_DETAIL_PARTITION_ELEMENT_INFO_H

// LOD includes.
#include <CGAL/Level_of_detail/Enumerations.h>

namespace CGAL {

	namespace Level_of_detail {

        namespace LOD = CGAL::Level_of_detail;

        class Partition_element_info {
        
        public:
            using Visibility_label = LOD::Visibility_label;

            Partition_element_info() :
            m_visibility_label(Visibility_label::OUTSIDE)
            { }

            inline Visibility_label& visibility_label() {
				return m_visibility_label;
			}

			inline const Visibility_label& visibility_label() const {
				return m_visibility_label;
			}

        private:
            Visibility_label m_visibility_label;
        };

    } // Level_of_detail

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_PARTITION_ELEMENT_INFO_H