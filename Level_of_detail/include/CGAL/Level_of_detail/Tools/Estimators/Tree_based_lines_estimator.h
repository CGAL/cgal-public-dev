#ifndef CGAL_LEVEL_OF_DETAIL_TREE_BASED_LINES_ESTIMATOR_H
#define CGAL_LEVEL_OF_DETAIL_TREE_BASED_LINES_ESTIMATOR_H

// STL includes.
#include <map>

// LOD includes.
#include <CGAL/Level_of_detail/Tools/Fitters/Line_to_points_fitter.h>

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

            using Neighbours       = typename Tree::Neighbours;
            using Point_identifier = typename Tree::Point_identifier;
            
            using Line_to_points_fitter = LOD::Line_to_points_fitter<Kernel>;
            using Elements_iterator     = typename Elements::const_iterator;

            using Lines_2 = std::map<Point_identifier, Line_2>;
            using Scores  = std::map<Point_identifier, FT>;

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
                for (Elements_iterator element = m_elements.begin(); element != m_elements.end(); ++element) {
                    
                    const Point_2 &point = get(m_point_map, *element);
                    m_tree.search_2(point, neighbours);

                     m_scores[*element] = line_to_points_fitter.fit_line_2(neighbours, m_tree.point_map(), line);
                    m_lines_2[*element] = line;
                }
            }
        };

    } // Level_of_detail

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_TREE_BASED_LINES_ESTIMATOR_H