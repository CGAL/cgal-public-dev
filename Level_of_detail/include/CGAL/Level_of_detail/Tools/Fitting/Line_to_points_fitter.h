#ifndef CGAL_LEVEL_OF_DETAIL_LINE_TO_POINTS_FITTER_H
#define CGAL_LEVEL_OF_DETAIL_LINE_TO_POINTS_FITTER_H

// STL includes.
#include <vector>

// CGAL includes.
#include <CGAL/Eigen_diagonalize_traits.h>
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
			using Local_FT 	    = typename Local_kernel::FT;
			using Local_point_2 = typename Local_kernel::Point_2;
			using Local_line_2  = typename Local_kernel::Line_2;

			using Diagonalize_traits = CGAL::Eigen_diagonalize_traits<Local_FT, 2>;

			template<class Elements, class Point_map>
			FT fit_line_2(const Elements &elements, const Point_map &point_map, Line_2 &line) const {
				
				CGAL_precondition(elements.size() > 0);
				using Const_elements_iterator = typename Elements::const_iterator;

				size_t i = 0;
				std::vector<Local_point_2> points(elements.size()); 

				Local_FT cx = Local_FT(0), cy = Local_FT(0);
				for (Const_elements_iterator ce_it = elements.begin(); ce_it != elements.end(); ++ce_it, ++i) {
					const Point_2 &point = get(point_map, *ce_it);

					const Local_FT x = static_cast<Local_FT>(CGAL::to_double(point.x()));
					const Local_FT y = static_cast<Local_FT>(CGAL::to_double(point.y()));

					points[i] = Local_point_2(x, y);

					cx += x;
					cy += y;
				}
				const Local_FT size = static_cast<Local_FT>(i);

				cx /= size;
				cy /= size;

				Local_point_2 centroid(cx, cy);
				Local_line_2  fitted_line;

				const FT quality = static_cast<FT>(CGAL::linear_least_squares_fitting_2(points.begin(), points.end(), fitted_line, centroid, CGAL::Dimension_tag<0>(), Local_kernel(), Diagonalize_traits()));
				line = Line_2(static_cast<FT>(fitted_line.a()), static_cast<FT>(fitted_line.b()), static_cast<FT>(fitted_line.c()));

				return quality;
			}
		};
	
	} // Level_of_detail

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_LINE_TO_POINTS_FITTER_H