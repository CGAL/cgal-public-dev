#ifndef CGAL_LEVEL_OF_DETAIL_LINEARITY_BASED_SORTING_2_H
#define CGAL_LEVEL_OF_DETAIL_LINEARITY_BASED_SORTING_2_H

// STL includes.
#include <map>
#include <list>

// CGAL includes.
#include <CGAL/Kd_tree.h>
#include <CGAL/Fuzzy_sphere.h>
#include <CGAL/Search_traits_2.h>

// LOD includes.
#include <CGAL/Level_of_detail/Tools/Fitters/Line_to_points_fitter.h>

namespace CGAL {

	namespace Level_of_detail {

        namespace LOD = CGAL::Level_of_detail;

		template<class PointIdentifier, class InputKernel, class InputElements, class PointMap>
		class Linearity_based_sorting_2 {

        public:
            using Point_identifier = PointIdentifier;
            using Kernel           = InputKernel;
            using Elements         = InputElements;
			using Point_map        = PointMap;

            using FT      = typename Kernel::FT;
            using Point_2 = typename Kernel::Point_2;
            using Line_2  = typename Kernel::Line_2;

            using Neighbours = std::list<Point_2>;
            using Scores     = std::map<Point_identifier, FT>;

			using Search_traits_2 = CGAL::Search_traits_2<Kernel>;
			using Search_circle   = CGAL::Fuzzy_sphere<Search_traits_2>;
			using Search_tree     = CGAL::Kd_tree<Search_traits_2>;

            using Identity_point_map    = CGAL::Identity_property_map<Point_2>;
            using Line_to_points_fitter = LOD::Line_to_points_fitter<Kernel>;

            Linearity_based_sorting_2(const Elements &elements, const Point_map &point_map, const FT local_search_radius) : 
            m_elements(elements),
            m_point_map(point_map),
            m_local_search_radius(local_search_radius) { 

                compute_scores();
            }

            bool operator() (const Point_identifier &i, const Point_identifier &j) const {
			    return m_scores.at(i) > m_scores.at(j);
			}

        private:
            const Elements  &m_elements;
            const Point_map &m_point_map;

            const FT m_local_search_radius;
            Scores   m_scores;

            void compute_scores() {
                m_scores.clear();

				Search_tree tree;
                create_tree(tree);

                for (typename Elements::const_iterator element = m_elements.begin(); element != m_elements.end(); ++element) {
                    const Point_2 &point = get(m_point_map, *element);

                    const FT score = compute_score(tree, point);
                    m_scores[*element] = score;
                }
            }

            void create_tree(Search_tree &tree) const {
                
                tree.clear();
                for (typename Elements::const_iterator element = m_elements.begin(); element != m_elements.end(); ++element) {
                    
                    const Point_2 &point = get(m_point_map, *element);
                    tree.insert(point);
                }
            }

            FT compute_score(const Search_tree &tree, const Point_2 &query) const {

                Neighbours neighbours;
				Search_circle circle(query, m_local_search_radius);
				tree.search(std::back_inserter(neighbours), circle);

                Line_2 line;
                Identity_point_map identity_point_map;

                Line_to_points_fitter line_to_points_fitter;
                return line_to_points_fitter.fit_line_2(neighbours, identity_point_map, line);
            }
        };

    } // Level_of_detail

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_LINEARITY_BASED_SORTING_2_H