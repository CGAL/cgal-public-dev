#ifndef CGAL_LEVEL_OF_DETAIL_PARTITION_ELEMENT_H
#define CGAL_LEVEL_OF_DETAIL_PARTITION_ELEMENT_H

// STL includes.
#include <vector>

namespace CGAL {

	namespace Level_of_detail {

        template<class InputKernel, class InputContainer>
        class Partition_element {
        
        public:
            using Kernel    = InputKernel;
            using Container = InputContainer;

            using FT = typename Kernel::FT;
            using Constraints = std::vector<bool>;

            using Const_iterator = typename Container::Vertex_const_iterator;
            
            template<class Elements, class Point_map>
            Partition_element(const Elements &elements, const Point_map &point_map) : 
            m_building_interior(FT(1)) {
                
                m_container.clear();
                using Const_elements_iterator = typename Elements::const_iterator;
                
                for (Const_elements_iterator ce_it = elements.begin(); ce_it != elements.end(); ++ce_it)
                    m_container.push_back(get(point_map, *ce_it));
            }

            template<class Point>
            inline bool locate(const Point &query) const {
                return m_container.has_on_bounded_side(query);
            }

            inline Constraints& constraints() {
				return m_constraints;
			}

			inline const Constraints& constraints() const {
				return m_constraints;
			}

            inline FT& building_interior() {
				return m_building_interior;
			}

			inline const FT& building_interior() const {
				return m_building_interior;
			}

            inline const Const_iterator begin() const {
                return m_container.vertices_begin();
            }

            inline const Const_iterator end() const {
                return m_container.vertices_end();
            }

            inline size_t size() const {
                return m_container.size();
            }

        private:
            Container   m_container;
            Constraints m_constraints;

            FT m_building_interior;
        };

    } // Level_of_detail

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_PARTITION_ELEMENT_H