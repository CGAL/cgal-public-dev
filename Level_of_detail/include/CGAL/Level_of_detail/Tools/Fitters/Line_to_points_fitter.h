#ifndef CGAL_LEVEL_OF_DETAIL_LINE_TO_POINTS_FITTER_H
#define CGAL_LEVEL_OF_DETAIL_LINE_TO_POINTS_FITTER_H

// STL includes.
#include <vector>

// CGAL includes.
#include <CGAL/linear_least_squares_fitting_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

namespace CGAL {

	namespace Level_of_detail {

		template<class InputKernel>
		class Line_to_points_fitter {

		public:
			using Kernel = InputKernel;
			
			using FT      = typename Kernel::FT;
			using Point_2 = typename Kernel::Point_2;
			using Line_2  = typename Kernel::Line_2;

			using Local_kernel  = CGAL::Exact_predicates_inexact_constructions_kernel;
			using Local_point_2 = typename Local_kernel::Point_2;
			using Local_line_2  = typename Local_kernel::Line_2;

			template<class Elements, class Point_map>
			FT fit_line_2(const Elements &elements, const Point_map &point_map, Line_2 &line) const {
				CGAL_precondition(elements.size() > 0);

				size_t i = 0;
				std::vector<Local_point_2> points(elements.size()); 

				for (typename Elements::const_iterator element = elements.begin(); element != elements.end(); ++element, ++i) {
					const Point_2 &point = get(point_map, *element);

					const double x = CGAL::to_double(point.x());
					const double y = CGAL::to_double(point.y());

					points[i] = Local_point_2(x, y);
				}

				Local_line_2 fitted_line;
				const FT quality = static_cast<FT>(CGAL::linear_least_squares_fitting_2(points.begin(), points.end(), fitted_line, CGAL::Dimension_tag<0>()));
				line = Line_2(static_cast<FT>(fitted_line.a()), static_cast<FT>(fitted_line.b()), static_cast<FT>(fitted_line.c()));

				return quality;
			}
		};
	
	} // Level_of_detail

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_LINE_TO_POINTS_FITTER_H