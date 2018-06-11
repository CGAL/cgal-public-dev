#ifndef CGAL_LEVEL_OF_DETAIL_PARTITION_ELEMENT_H
#define CGAL_LEVEL_OF_DETAIL_PARTITION_ELEMENT_H

// LOD includes.
#include <CGAL/Level_of_detail/Partitioning/Data_structures/Partition_element_info.h>

namespace CGAL {

	namespace Level_of_detail {

        namespace LOD = CGAL::Level_of_detail;

        template<class InputKernel, class InputContainer>
        class Partition_element {
        
        public:
            using Kernel    = InputKernel;
            using Container = InputContainer;

            using FT = typename Kernel::FT;
            using const_iterator = typename Container::Vertex_const_iterator;

            using Partition_element_info = LOD::Partition_element_info;
            
            template<class Elements, class Point_map>
            Partition_element(const Elements &elements, const Point_map &point_map) {
                
                m_container.clear();
                using Const_elements_iterator = typename Elements::const_iterator;
                
                for (Const_elements_iterator ce_it = elements.begin(); ce_it != elements.end(); ++ce_it)
                    m_container.push_back(get(point_map, *ce_it));
            }

            template<class Point>
            inline bool has_on_bounded_side(const Point &query) const {
                return m_container.has_on_bounded_side(query);
            }

            inline const const_iterator begin() const {
                return m_container.vertices_begin();
            }

            inline const const_iterator end() const {
                return m_container.vertices_end();
            }

            inline Partition_element_info& info() {
				return m_partition_element_info;
			}

			inline const Partition_element_info& info() const {
				return m_partition_element_info;
			}

            inline size_t size() const {
                return m_container.size();
            }

        private:
            Container              m_container;
            Partition_element_info m_partition_element_info;
        };

    } // Level_of_detail

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_PARTITION_ELEMENT_H