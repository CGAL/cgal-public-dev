#ifndef CGAL_LEVEL_OF_DETAIL_PLANE_TO_POINTS_FITTER_H
#define CGAL_LEVEL_OF_DETAIL_PLANE_TO_POINTS_FITTER_H

// STL includes.
#include <vector>

// CGAL includes.
#include <CGAL/Eigen_diagonalize_traits.h>
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
			using Local_FT 		= typename Local_kernel::FT;
			using Local_point_3 = typename Local_kernel::Point_3;
			using Local_plane_3 = typename Local_kernel::Plane_3;

			using Diagonalize_traits = CGAL::Eigen_diagonalize_traits<Local_FT, 3>;

			template<class Elements, class Point_map>
			FT fit_plane(const Elements &elements, const Point_map &point_map, Plane_3 &plane) const {
				
				CGAL_precondition(elements.size() > 2);
				using Const_elements_iterator = typename Elements::const_iterator;

				size_t i = 0;
				std::vector<Local_point_3> points(elements.size());
				
				Local_FT cx = Local_FT(0), cy = Local_FT(0), cz = Local_FT(0);
				for (Const_elements_iterator element = elements.begin(); element != elements.end(); ++element, ++i) {
					const Point_3 &point = get(point_map, *element);

					const Local_FT x = static_cast<Local_FT>(CGAL::to_double(point.x()));
					const Local_FT y = static_cast<Local_FT>(CGAL::to_double(point.y()));
					const Local_FT z = static_cast<Local_FT>(CGAL::to_double(point.z()));

					points[i] = Local_point_3(x, y, z);

					cx += x;
					cy += y;
					cz += z;
				}
				const Local_FT size = static_cast<Local_FT>(i);

				cx /= size;
				cy /= size;
				cz /= size;

				Local_point_3 centroid(cx, cy, cz);
				Local_plane_3 fitted_plane;

				const FT quality = static_cast<FT>(CGAL::linear_least_squares_fitting_3(points.begin(), points.end(), fitted_plane, centroid, CGAL::Dimension_tag<0>(), Local_kernel(), Diagonalize_traits()));
				plane = Plane_3(static_cast<FT>(fitted_plane.a()), static_cast<FT>(fitted_plane.b()), static_cast<FT>(fitted_plane.c()), static_cast<FT>(fitted_plane.d()));

				return quality;
			}
		};
	
	} // Level_of_detail

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_PLANE_TO_POINTS_FITTER_H