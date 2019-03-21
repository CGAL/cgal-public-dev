// Copyright (c) 2019 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is a part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY, AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
//
// Author(s)     : Dmitry Anisimov, Simon Giraudot, Pierre Alliez, Florent Lafarge, and Andreas Fabri

#ifndef CGAL_LEVELS_OF_DETAIL_INTERNAL_UTILS_H
#define CGAL_LEVELS_OF_DETAIL_INTERNAL_UTILS_H

#include <CGAL/license/Levels_of_detail.h>

// STL includes.
#include <map>
#include <vector>

// CGAL includes.
#include <CGAL/assertions.h>
#include <CGAL/linear_least_squares_fitting_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

// Internal includes.
#include <CGAL/Levels_of_detail/internal/number_utils.h>

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

  template<typename Point>
  class Indexer {
  
  public:
    std::size_t operator()(const Point& point) {
      const auto pair = m_indices.insert(
        std::make_pair(
          point, 
          m_indices.size()));
      const auto& item = pair.first;
      const std::size_t idx = item->second;
      return idx;
    }
    void clear() { m_indices.clear(); }

  private:
    std::map<Point, std::size_t> m_indices;
  };

  template<
  typename Point_2, 
  typename Plane_3>
  typename Kernel_traits<Point_2>::Kernel::Point_3
  position_on_plane_3(
    const Point_2& point, 
    const Plane_3& plane) {

    using Traits = typename Kernel_traits<Point_2>::Kernel;
    using FT = typename Traits::FT;
    using Line_3 = typename Traits::Line_3;
    using Point_3 = typename Traits::Point_3;
    using Vector_3 = typename Traits::Vector_3;
    using Intersect_3 = typename Traits::Intersect_3;

    static Vector_3 vertical(FT(0), FT(0), FT(1));
    const Line_3 line(Point_3(point.x(), point.y(), FT(0)), vertical);
    typename CGAL::cpp11::result_of<Intersect_3(Line_3, Plane_3)>::type
      inter = CGAL::intersection(line, plane);
    if (inter)
      if (const Point_3* p = boost::get<Point_3>(&*inter))
        return *p;
    
    std::cerr << 
      "Error (position_on_plane): cannot compute the 3D position!" 
    << std::endl;
    return Point_3(FT(0), FT(0), FT(0));
  }

  template<
  typename Item_range, 
  typename Point_map_3, 
  typename Plane_3>
	typename Kernel_traits<
  typename boost::property_traits<Point_map_3>::value_type>::Kernel::FT
  plane_from_points_3(
    const Item_range& item_range, 
    const Point_map_3 point_map_3, 
    Plane_3& plane) {

    using Traits = typename Kernel_traits<Plane_3>::Kernel;
    using FT = typename Traits::FT;
    using Point_3 = typename Traits::Point_3;

    using Local_traits 
    = CGAL::Exact_predicates_inexact_constructions_kernel;
		using Local_FT = typename Local_traits::FT;
		using Local_point_3 = typename Local_traits::Point_3;
		using Local_plane_3 = typename Local_traits::Plane_3;

		CGAL_assertion(item_range.size() > 0);
		std::vector<Local_point_3> points;
    points.reserve(item_range.size());
				
		for (std::size_t i = 0; i < item_range.size(); ++i) {
			const Point_3& p = get(point_map_3, *(item_range.begin() + i));

			const Local_FT x = static_cast<Local_FT>(CGAL::to_double(p.x()));
			const Local_FT y = static_cast<Local_FT>(CGAL::to_double(p.y()));
			const Local_FT z = static_cast<Local_FT>(CGAL::to_double(p.z()));

			points.push_back(Local_point_3(x, y, z));
		}
    CGAL_assertion(points.size() == item_range.size());

		Local_plane_3 fitted_plane;
    Local_point_3 fitted_centroid;

		const FT quality = static_cast<FT>(
      CGAL::linear_least_squares_fitting_3(
        points.begin(), points.end(), 
        fitted_plane, fitted_centroid, 
        CGAL::Dimension_tag<0>()));

		plane = Plane_3(
      static_cast<FT>(fitted_plane.a()), 
      static_cast<FT>(fitted_plane.b()), 
      static_cast<FT>(fitted_plane.c()), 
      static_cast<FT>(fitted_plane.d()));

    return quality;
	}

  template<
  typename Item_range, 
  typename Point_map_2,
  typename Point_2>
  void bounding_box_2(
    const Item_range& item_range,
    const Point_map_2 point_map_2,
    std::vector<Point_2>& bbox) {
    
    using Traits = typename Kernel_traits<Point_2>::Kernel;
    using FT = typename Traits::FT;

    CGAL_assertion(item_range.size() > 0);
    FT minx = internal::max_value<FT>(), miny = internal::max_value<FT>();
    FT maxx = -internal::max_value<FT>(), maxy = -internal::max_value<FT>();

    for (std::size_t i = 0; i < item_range.size(); ++i) {
      const Point_2& p = get(point_map_2, *(item_range.begin() + i));
      minx = CGAL::min(minx, p.x()); miny = CGAL::min(miny, p.y());
      maxx = CGAL::max(maxx, p.x()); maxy = CGAL::max(maxy, p.y());
    }

    bbox.clear(); bbox.reserve(4);
    bbox.push_back(Point_2(minx, miny)); bbox.push_back(Point_2(maxx, miny));
    bbox.push_back(Point_2(maxx, maxy)); bbox.push_back(Point_2(minx, maxy));
    CGAL_assertion(bbox.size() == 4);
  }

  template<typename Point_3>
  typename Kernel_traits<Point_3>::Kernel::Point_2
  point_2_from_point_3(const Point_3& point_3) {
    return typename Kernel_traits<Point_3>::Kernel::Point_2(
      point_3.x(), point_3.y());
  }

  template<
  typename Face_handle,
  typename Triangle_3>
  void triangle_3(
    const Face_handle& fh, 
    Triangle_3& triangle) {
    
    using Traits = typename Kernel_traits<Triangle_3>::Kernel;
    using Point_3 = typename Traits::Point_3;
    
    Point_3 p0, p1, p2;
    point_3(fh, 0, p0);
    point_3(fh, 1, p1);
    point_3(fh, 2, p2);
    triangle = Triangle_3(p0, p1, p2);
  }

  template<
  typename Face_handle,
  typename Point_3>
  void point_3(
    const Face_handle& fh, 
    const std::size_t idx,
    Point_3& p) {
    p = Point_3(fh->vertex(idx)->point().x(), 
                fh->vertex(idx)->point().y(), 
                fh->info().z[idx]);
  }

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_INTERNAL_UTILS_H
