#ifndef CGAL_LEVEL_OF_DETAIL_KD_TREE_WITH_DATA_CREATOR_H
#define CGAL_LEVEL_OF_DETAIL_KD_TREE_WITH_DATA_CREATOR_H

// STL includes.
#include <list>
#include <utility>

// CGAL includes.
#include <CGAL/Kd_tree.h>
#include <CGAL/property_map.h>
#include <CGAL/Fuzzy_sphere.h>
#include <CGAL/Search_traits_2.h>
#include <CGAL/Search_traits_adapter.h>

namespace CGAL {

	namespace Level_of_detail {

        template<class InputKernel, class ElementIdentifier, class InputElements, class PointMap>
		class Kd_tree_with_data_creator {

        public:
            using Kernel             = InputKernel;
            using Element_identifier = ElementIdentifier;
            using Elements           = InputElements;
			using Point_map          = PointMap;

            using FT      = typename Kernel::FT;
            using Point_2 = typename Kernel::Point_2;

            using Search_traits_2 = CGAL::Search_traits_2<Kernel>;
            
            using Tree_element = std::pair<Point_2, Element_identifier>;
            using Neighbours   = std::list<Tree_element>;

            using Tree_point_map            = CGAL::First_of_pair_property_map<Tree_element>;
            using Tree_point_identifier_map = CGAL::Second_of_pair_property_map<Tree_element>;
            
            using Search_traits           = CGAL::Search_traits_adapter<Tree_element, Tree_point_map, Search_traits_2>;
            using Const_elements_iterator = typename Elements::const_iterator;
            
            /*
            using Tree_element = Point_identifier;
            using Neighbours   = std::list<Tree_element>;
            
            using Tree_point_map            = Point_map;
            using Tree_point_identifier_map = CGAL::Identity_property_map<Point_identifier>;
            
            using Search_traits = CGAL::Search_traits_adapter<Tree_element, Tree_point_map, Search_traits_2>;
            using Splitter      = typename Search_tree::Splitter; */

			using Search_circle   = CGAL::Fuzzy_sphere<Search_traits>;
			using Search_tree     = CGAL::Kd_tree<Search_traits>;

            Kd_tree_with_data_creator(const Elements &elements, const Point_map &point_map, const FT local_search_radius) : 
            m_elements(elements),
            m_point_map(point_map),
            /* m_search_traits(m_point_map), */
            m_local_search_radius(local_search_radius) { 

                create_tree_2();
            }

            void set_local_search_radius(const FT local_search_radius) {
                
                CGAL_precondition(local_search_radius > FT(0));
                m_local_search_radius = local_search_radius;
            }

            void search_2(const Point_2 &query, Neighbours &neighbours) const {
                neighbours.clear();

                /* Search_circle circle(query, m_local_search_radius, FT(0), m_search_traits); */

                Search_circle circle(query, m_local_search_radius);
				m_tree.search(std::back_inserter(neighbours), circle);
            }

            inline const Tree_point_map& point_map() const {
                return m_tree_point_map;
            }

            inline const Tree_point_identifier_map& point_identifier_map() const {
                return m_tree_point_identifier_map;
            }

        private:
            const Elements  &m_elements;
            const Point_map &m_point_map;

            Search_tree               m_tree;
            Tree_point_map            m_tree_point_map;
            Tree_point_identifier_map m_tree_point_identifier_map;

            /*
            const Splitter      m_splitter;
            const Search_traits m_search_traits; */

            FT m_local_search_radius;

            void create_tree_2() {
                /* m_tree = Search_tree(m_elements.begin(), m_elements.end(), m_splitter, m_search_traits); */
                
                m_tree.clear();
                for (Const_elements_iterator ce_it = m_elements.begin(); ce_it != m_elements.end(); ++ce_it) {
                    
                    const Point_2 &point = get(m_point_map, *ce_it);
                    m_tree.insert(std::make_pair(point, *ce_it));
                }
            }
        };

    } // Level_of_detail

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_KD_TREE_WITH_DATA_CREATOR_H