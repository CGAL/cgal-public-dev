#ifndef CGAL_LEVELS_OF_DETAIL_UTILITIES_H
#define CGAL_LEVELS_OF_DETAIL_UTILITIES_H

// STL includes.
#include <vector>

// CGAL includes.
#include <CGAL/Eigen_diagonalize_traits.h>
#include <CGAL/linear_least_squares_fitting_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

namespace CGAL {

namespace Levels_of_detail {

namespace internal {

	template<class Items, class Point_map, class Plane_3>
	void fit_plane_to_points_3(
    const Items &items, 
    const Point_map point_map, 
    Plane_3 &plane) {

    using Traits = 
    typename Kernel_traits<
    typename boost::property_traits<Point_map>::value_type>::Kernel;
    
    using FT = typename Traits::FT;
    using Point_3 = typename Traits::Point_3;

    using Local_traits 
    = CGAL::Exact_predicates_inexact_constructions_kernel;

		using Local_FT = typename Local_traits::FT;
		using Local_point_3 = typename Local_traits::Point_3;
		using Local_plane_3 = typename Local_traits::Plane_3;

		using Diagonalize_traits = CGAL::Eigen_diagonalize_traits<Local_FT, 3>;

		CGAL_precondition(items.size() > 2);
		std::vector<Local_point_3> points(items.size());
				
    std::size_t i = 0;
		for (auto it = items.begin(); it != items.end(); ++it, ++i) {
			const Point_3 &point = get(point_map, *it);

			const Local_FT x = static_cast<Local_FT>(CGAL::to_double(point.x()));
			const Local_FT y = static_cast<Local_FT>(CGAL::to_double(point.y()));
			const Local_FT z = static_cast<Local_FT>(CGAL::to_double(point.z()));

			points[i] = Local_point_3(x, y, z);
		}

		Local_plane_3 fitted_plane;
    Local_point_3 fitted_centroid;

		CGAL::linear_least_squares_fitting_3(
      points.begin(), points.end(), 
      fitted_plane, fitted_centroid, CGAL::Dimension_tag<0>(),
      Local_traits(), Diagonalize_traits());

		plane = Plane_3(
      static_cast<FT>(fitted_plane.a()), 
      static_cast<FT>(fitted_plane.b()), 
      static_cast<FT>(fitted_plane.c()), 
      static_cast<FT>(fitted_plane.d()));
	}

  template<class Items, class Point_map, class Plane_3, class Point_3>
  void compute_bounding_box_3(
    const Items &items, 
    const Point_map &point_map,
    const Plane_3 &plane, 
    std::vector<Point_3> &bounding_box) {

    using Traits = 
    typename Kernel_traits<
    typename boost::property_traits<Point_map>::value_type>::Kernel;
    using FT = typename Traits::FT;
                
    CGAL_precondition(items.size() > 0);

    FT minx =  std::numeric_limits<FT>::max(); 
    FT miny =  std::numeric_limits<FT>::max();
    FT maxx = -std::numeric_limits<FT>::max();
    FT maxy = -std::numeric_limits<FT>::max();
                
    FT z = FT(0), size = FT(0);
    for (auto it = items.begin(); it != items.end(); ++it, size += FT(1)) {
					
      const Point_3 &point = get(point_map, *it);
      const Point_3 projected = plane.projection(point);

      minx = CGAL::min(minx, projected.x());
      miny = CGAL::min(miny, projected.y());

      maxx = CGAL::max(maxx, projected.x());
      maxy = CGAL::max(maxy, projected.y());

      z += projected.z();
    }
    z /= size;

    bounding_box.clear();
    bounding_box.resize(4);

    bounding_box[0] = Point_3(minx, miny, z);
    bounding_box[1] = Point_3(maxx, miny, z);
    bounding_box[2] = Point_3(maxx, maxy, z);
    bounding_box[3] = Point_3(minx, maxy, z);
  }

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_UTILITIES_H
