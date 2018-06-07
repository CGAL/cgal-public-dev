#ifndef CGAL_LEVEL_OF_DETAIL_POINTS_TO_LINE_PROJECTOR_H
#define CGAL_LEVEL_OF_DETAIL_POINTS_TO_LINE_PROJECTOR_H

namespace CGAL {

	namespace Level_of_detail {

		template<class InputKernel>
		class Points_to_line_projector {

		public:
			using Kernel = InputKernel;
            
            using Line_2  = typename Kernel::Line_2;
            using Point_2 = typename Kernel::Point_2;

            template<class Elements, class Point_map, class Output>
            void project_on_line_2(const Elements &elements, const Point_map &point_map, const Line_2 &line, Output &output) const {
                
                using Const_elements_iterator = typename Elements::const_iterator;
                output.clear();

                for (Const_elements_iterator element = elements.begin(); element != elements.end(); ++element) {	
                    const Point_2 &point = get(point_map, *element);

                    output.push_back(line.projection(point));
                }
            }
        };

    } // Level_of_detail

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_POINTS_TO_LINE_PROJECTOR_H