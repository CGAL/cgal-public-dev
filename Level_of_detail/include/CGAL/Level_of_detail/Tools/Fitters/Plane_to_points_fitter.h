#ifndef CGAL_LEVEL_OF_DETAIL_PLANE_TO_POINTS_FITTER_H
#define CGAL_LEVEL_OF_DETAIL_PLANE_TO_POINTS_FITTER_H

// STL includes.
#include <vector>

// CGAL includes.
#include <CGAL/linear_least_squares_fitting_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

namespace CGAL {

	namespace Level_of_detail {

		template<class InputKernel>
		class Plane_to_points_fitter {

		public:
			using Kernel = InputKernel;
			
			using FT      = typename Kernel::FT;
			using Point_3 = typename Kernel::Point_3;
			using Plane_3 = typename Kernel::Plane_3;

			using Local_kernel  = CGAL::Exact_predicates_inexact_constructions_kernel;
			using Local_point_3 = typename Local_kernel::Point_3;
			using Local_plane_3 = typename Local_kernel::Plane_3;

			template<class Elements, class Point_map>
			FT fit_plane(const Elements &elements, const Point_map &point_map, Plane_3 &plane) const {
				CGAL_precondition(elements.size() > 2);

				size_t i = 0;
				std::vector<Local_point_3> points(elements.size()); 

				for (typename Elements::const_iterator element = elements.begin(); element != elements.end(); ++element, ++i) {
					const Point_3 &point = get(point_map, *element);

					const double x = CGAL::to_double(point.x());
					const double y = CGAL::to_double(point.y());
					const double z = CGAL::to_double(point.z());

					points[i] = Local_point_3(x, y, z);
				}

				Local_plane_3 fitted_plane;
				const FT quality = static_cast<FT>(CGAL::linear_least_squares_fitting_3(points.begin(), points.end(), fitted_plane, CGAL::Dimension_tag<0>()));
				plane = Plane_3(static_cast<FT>(fitted_plane.a()), static_cast<FT>(fitted_plane.b()), static_cast<FT>(fitted_plane.c()), static_cast<FT>(fitted_plane.d()));

				return quality;
			}
		};
	
	} // Level_of_detail

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_PLANE_TO_POINTS_FITTER_H