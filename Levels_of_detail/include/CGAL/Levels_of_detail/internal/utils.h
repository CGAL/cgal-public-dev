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
//
// Author(s)     : Dmitry Anisimov, Simon Giraudot, Pierre Alliez, Florent Lafarge, and Andreas Fabri
//

#ifndef CGAL_LEVELS_OF_DETAIL_INTERNAL_UTILS_H
#define CGAL_LEVELS_OF_DETAIL_INTERNAL_UTILS_H

#include <CGAL/license/Levels_of_detail.h>

// STL includes.
#include <map>
#include <vector>

// CGAL includes.
#include <CGAL/utils.h>
#include <CGAL/assertions.h>
#include <CGAL/Kernel/global_functions.h>
#include <CGAL/compute_average_spacing.h>
#include <CGAL/linear_least_squares_fitting_2.h>
#include <CGAL/linear_least_squares_fitting_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Barycentric_coordinates_2/Segment_coordinates_2.h>

// Internal includes.
#include <CGAL/Levels_of_detail/internal/number_utils.h>

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

  template<typename FT>
  struct Compare_scores {
    const std::vector<FT>& m_scores;
      
    Compare_scores(const std::vector<FT>& scores) : 
    m_scores(scores) 
    { }

    bool operator()(const std::size_t i, const std::size_t j) const {
      CGAL_precondition(i >= 0 && i < m_scores.size());
      CGAL_precondition(j >= 0 && j < m_scores.size());
      return m_scores[i] > m_scores[j];
    }
  };

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
  typename Point_map_2, 
  typename Line_2>
	typename Kernel_traits<
  typename boost::property_traits<Point_map_2>::value_type>::Kernel::FT
  line_from_points_2(
    const Item_range& item_range, 
    const Point_map_2& point_map_2, 
    Line_2& line) {

    using Traits = typename Kernel_traits<Line_2>::Kernel;
    using FT = typename Traits::FT;
    using Point_2 = typename Traits::Point_2;

    using Local_traits 
    = CGAL::Exact_predicates_inexact_constructions_kernel;
		using Local_FT = typename Local_traits::FT;
		using Local_line_2 = typename Local_traits::Line_2;
    using Local_point_2 = typename Local_traits::Point_2;

		CGAL_assertion(item_range.size() > 0);
		std::vector<Local_point_2> points;
    points.reserve(item_range.size());
				
		for (std::size_t i = 0; i < item_range.size(); ++i) {
			const Point_2& p = get(point_map_2, *(item_range.begin() + i));

			const Local_FT x = static_cast<Local_FT>(CGAL::to_double(p.x()));
			const Local_FT y = static_cast<Local_FT>(CGAL::to_double(p.y()));

			points.push_back(Local_point_2(x, y));
		}
    CGAL_assertion(points.size() == item_range.size());

		Local_line_2 fitted_line;
    Local_point_2 fitted_centroid;

		const FT quality = static_cast<FT>(
      CGAL::linear_least_squares_fitting_2(
        points.begin(), points.end(), 
        fitted_line, fitted_centroid, 
        CGAL::Dimension_tag<0>()));

		line = Line_2(
      static_cast<FT>(fitted_line.a()), 
      static_cast<FT>(fitted_line.b()), 
      static_cast<FT>(fitted_line.c()));

    return quality;
  }

  template<
  typename Item_range, 
  typename Point_map_3, 
  typename Plane_3>
	typename Kernel_traits<
  typename boost::property_traits<Point_map_3>::value_type>::Kernel::FT
  plane_from_points_3(
    const Item_range& item_range, 
    const Point_map_3& point_map_3, 
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
    const Point_map_2& point_map_2,
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
  typename Segment_2, 
  typename Point_2>
  void bounding_box_2(
    const std::vector<Segment_2>& segments,
    std::vector<Point_2>& bbox) {
    
    using Traits = typename Kernel_traits<Segment_2>::Kernel;
    using FT = typename Traits::FT;

    CGAL_assertion(segments.size() > 0);
    FT minx = internal::max_value<FT>(), miny = internal::max_value<FT>();
    FT maxx = -internal::max_value<FT>(), maxy = -internal::max_value<FT>();

    for (const auto& segment : segments) {                      
      const Point_2& source = segment.source();
      const Point_2& target = segment.target();

      minx = CGAL::min(minx, source.x()); minx = CGAL::min(minx, target.x());
      miny = CGAL::min(miny, source.y()); miny = CGAL::min(miny, target.y());
      maxx = CGAL::max(maxx, source.x()); maxx = CGAL::max(maxx, target.x());
      maxy = CGAL::max(maxy, source.y()); maxy = CGAL::max(maxy, target.y());
    }

    bbox.clear(); bbox.reserve(4);
    bbox.push_back(Point_2(minx, miny)); bbox.push_back(Point_2(maxx, miny));
    bbox.push_back(Point_2(maxx, maxy)); bbox.push_back(Point_2(minx, maxy));
    CGAL_assertion(bbox.size() == 4);
  }

  template<typename Segment_2>
  typename Kernel_traits<Segment_2>::Kernel::FT
  average_spacing_2(
    const std::vector<Segment_2>& segments,
    const std::size_t num_neighbors) {

    using Traits = typename Kernel_traits<Segment_2>::Kernel;
    using FT = typename Traits::FT;
    using Point_2 = typename Traits::Point_2;

    using Local_traits = CGAL::Exact_predicates_inexact_constructions_kernel;
		using Local_FT = typename Local_traits::FT;
    using Local_point_3 = typename Local_traits::Point_3;

    CGAL_assertion(segments.size() > 0);
    std::vector<Local_point_3> points;
    points.reserve(segments.size() * 2);

    for (const auto& segment : segments) {
      const Point_2& source = segment.source();
      const Point_2& target = segment.target();

      const Local_FT sx = static_cast<Local_FT>(CGAL::to_double(source.x()));
      const Local_FT sy = static_cast<Local_FT>(CGAL::to_double(source.y()));
      const Local_FT tx = static_cast<Local_FT>(CGAL::to_double(target.x()));
      const Local_FT ty = static_cast<Local_FT>(CGAL::to_double(target.y()));

      points.push_back(Local_point_3(sx, sy, Local_FT(0)));
      points.push_back(Local_point_3(tx, ty, Local_FT(0)));
    }

    const Local_FT average_spacing = 
    CGAL::compute_average_spacing<CGAL::Sequential_tag>(
      points, num_neighbors, CGAL::parameters::point_map(
        CGAL::Identity_property_map<Local_point_3>()).
        geom_traits(Local_traits()));
                
    return static_cast<FT>(average_spacing);
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

  template<
  typename Point_2,
  typename Triangle_2,
  typename FT>
  bool is_within_triangle(
    const Point_2& query,
    const Triangle_2& triangle,
    const FT tol) {

    using Traits = typename Kernel_traits<Point_2>::Kernel;
    using Line_2 = typename Traits::Line_2;

    if (triangle.has_on_bounded_side(query) || 
        triangle.has_on_boundary(query)) 
      return true;
                
    for (std::size_t i = 0; i < 3; ++i) {
      const std::size_t ip = (i + 1) % 3;

      const Point_2& p1 = triangle.vertex(i);
      const Point_2& p2 = triangle.vertex(ip);
      const Line_2 line = Line_2(p1, p2);

      const Point_2 projected = line.projection(query);
      const FT squared_distance = CGAL::squared_distance(query, projected);

      const Traits traits;
      const auto pair = Barycentric_coordinates::
        compute_segment_coordinates_2(p1, p2, projected, traits);

      const FT squared_tolerance = tol * tol;
      const FT epst = FT(6) / FT(5);
      const FT epsb = -FT(1) / FT(5);

      if (pair[0] > epsb && pair[1] > epsb && 
          pair[0] < epst && pair[1] < epst && 
          squared_distance < squared_tolerance) 
        return true;
    }
    return false;
  }

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_INTERNAL_UTILS_H
