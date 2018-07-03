#ifndef CGAL_LEVEL_OF_DETAIL_TREE_BASED_LINES_ESTIMATOR_H
#define CGAL_LEVEL_OF_DETAIL_TREE_BASED_LINES_ESTIMATOR_H

// STL includes.
#include <map>

// LOD includes.
#include <CGAL/Level_of_detail/internal/Fitting/Line_to_points_fitter.h>

namespace CGAL {

	namespace Level_of_detail {

        namespace LOD = CGAL::Level_of_detail;

		template<class InputKernel, class InputElements, class PointMap, class TreeWrapper>
        class Tree_based_lines_estimator {

        public:
            using Kernel    = InputKernel;
            using Elements  = InputElements;
			using Point_map = PointMap;
            using Tree      = TreeWrapper;

            using FT      = typename Kernel::FT;
            using Point_2 = typename Kernel::Point_2;
            using Line_2  = typename Kernel::Line_2;

            using Neighbours         = typename Tree::Neighbours;
            using Element_identifier = typename Tree::Element_identifier;
            
            using Line_to_points_fitter   = LOD::Line_to_points_fitter<Kernel>;
            using Const_elements_iterator = typename Elements::const_iterator;

            using Lines_2 = std::map<Element_identifier, Line_2>;
            using Scores  = std::map<Element_identifier, FT>;

            Tree_based_lines_estimator(const Elements &elements, const Point_map &point_map, const Tree &tree) : 
            m_elements(elements),
            m_point_map(point_map), 
            m_tree(tree) { 

                estimate_lines_2();
            }

            inline const Scores& scores() const {
                return m_scores;
            }

            inline const Lines_2& lines_2() const {
                return m_lines_2;
            }


            class Sorter {

            public:
              Sorter (const Scores& scores) :
                m_scores(scores) 
              { }

              bool operator() (const Element_identifier& i, const Element_identifier& j) const
              {
                return m_scores.at(i) > m_scores.at(j);
              }

            private:
              const Scores& m_scores;
            };

            Sorter sorter() const
            {
              return Sorter (m_scores);
            }


        private:
            const Elements  &m_elements;
            const Point_map &m_point_map;
            const Tree      &m_tree;

            Scores  m_scores;
            Lines_2 m_lines_2;

            void estimate_lines_2() {
                
                m_scores.clear();
                m_lines_2.clear();
                
                Line_2     line;
                Neighbours neighbours;

                const Line_to_points_fitter line_to_points_fitter;
                for (Const_elements_iterator ce_it = m_elements.begin(); ce_it != m_elements.end(); ++ce_it) {
                    
                    const Point_2 &point = get(m_point_map, *ce_it);
                    m_tree.search_2(point, neighbours);

                     m_scores[*ce_it] = line_to_points_fitter.fit_line_2(neighbours, m_tree.point_map(), line);
                    m_lines_2[*ce_it] = line;
                }
            }
        };

    } // Level_of_detail

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_TREE_BASED_LINES_ESTIMATOR_H
