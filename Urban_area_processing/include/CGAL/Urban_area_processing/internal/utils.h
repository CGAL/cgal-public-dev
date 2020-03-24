// Copyright (c) 2020 SARL GeometryFactory (France).
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

#ifndef CGAL_URBAN_AREA_PROCESSING_INTERNAL_UTILS_H
#define CGAL_URBAN_AREA_PROCESSING_INTERNAL_UTILS_H

// #include <CGAL/license/Urban_area_processing.h>

// STL includes.
#include <vector>

// TODO:
// 1. Should I make some of these functions public like the one: boundary_points_on_line_2 e.g.?

namespace CGAL {
namespace Urban_area_processing {
namespace internal {

  template<typename GeomTraits> 
  class Default_sqrt {
    
  private:
    using Traits = GeomTraits;
    using FT = typename Traits::FT;

  public:
    FT operator()(const FT value) const { 
      
      CGAL_precondition(value >= FT(0));
      return static_cast<FT>(CGAL::sqrt(CGAL::to_double(value)));
    }
  };

  BOOST_MPL_HAS_XXX_TRAIT_NAMED_DEF(Has_nested_type_Sqrt, Sqrt, false)

  // Case: do_not_use_default = false.
  template<typename GeomTraits, 
  bool do_not_use_default = Has_nested_type_Sqrt<GeomTraits>::value>
  class Get_sqrt {
        
  public:
    using Traits = GeomTraits;
    using Sqrt = Default_sqrt<Traits>;

    static Sqrt sqrt_object(const Traits& ) { 
      return Sqrt();
    }
  };

  // Case: do_not_use_default = true.
  template<typename GeomTraits>
  class Get_sqrt<GeomTraits, true> {
        
  public:
    using Traits = GeomTraits;
    using Sqrt = typename Traits::Sqrt;

    static Sqrt sqrt_object(const Traits& traits) { 
      return traits.sqrt_object();
    }
  };

  template<typename FT>
  static FT max_value() {
    return FT(1000000000000);
  }

  template<typename FT>
  static FT tolerance() {
    return FT(1) / FT(100000);
  }

  template<
  typename Face_handle,
  typename Point_3>
  void point_3(
    const Face_handle& fh, 
    const std::size_t idx,
    Point_3& p) {
    p = Point_3(
      fh->vertex(idx)->point().x(), 
      fh->vertex(idx)->point().y(), 
      fh->info().z[idx]);
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

  template<typename Point_3>
  typename Kernel_traits<Point_3>::Kernel::Point_2
  point_2_from_point_3(const Point_3& point_3) {
    return typename Kernel_traits<Point_3>::Kernel::Point_2(
      point_3.x(), point_3.y());
  }

  template<typename Vector_3>
  typename Kernel_traits<Vector_3>::Kernel::FT
  angle_3d(
    const Vector_3& v1, 
    const Vector_3& v2) {
        
    using Traits = typename Kernel_traits<Vector_3>::Kernel;
    using FT = typename Traits::FT;

    const double a = CGAL::to_double(v1 * v2) / (
      CGAL::sqrt(CGAL::to_double(v1.squared_length())) * 
      CGAL::sqrt(CGAL::to_double(v2.squared_length())));

    if (a < -1.0) return static_cast<FT>(std::acos(-1.0) / CGAL_PI * 180.0);
    else if (a > 1.0) return static_cast<FT>(std::acos(1.0) / CGAL_PI * 180.0);
    return static_cast<FT>(std::acos(a) / CGAL_PI * 180.0);
  }

  template<
  typename Item_range, 
  typename Point_map_3, 
  typename Plane_3>
  typename Kernel_traits<Plane_3>::Kernel::FT
  plane_from_points_3(
    const Item_range& item_range, 
    const Point_map_3& point_map_3, 
    const std::vector<std::size_t>& indices,
    Plane_3& plane) {

    using Traits = typename Kernel_traits<Plane_3>::Kernel;
    using FT = typename Traits::FT;
    using Point_3 = typename Traits::Point_3;

    using Local_traits 
    = CGAL::Exact_predicates_inexact_constructions_kernel;
		using Local_FT = typename Local_traits::FT;
		using Local_point_3 = typename Local_traits::Point_3;
		using Local_plane_3 = typename Local_traits::Plane_3;

		CGAL_assertion(indices.size() > 0);
		std::vector<Local_point_3> points;
    points.reserve(indices.size());
				
		for (const std::size_t idx : indices) {
			const Point_3& p = get(point_map_3, *(item_range.begin() + idx));

			const Local_FT x = static_cast<Local_FT>(CGAL::to_double(p.x()));
			const Local_FT y = static_cast<Local_FT>(CGAL::to_double(p.y()));
			const Local_FT z = static_cast<Local_FT>(CGAL::to_double(p.z()));

			points.push_back(Local_point_3(x, y, z));
		}
    CGAL_assertion(points.size() == indices.size());

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
  typename Point_map_3, 
  typename Plane_3>
  typename Kernel_traits<Plane_3>::Kernel::FT
  plane_from_points_3(
    const Item_range& item_range, 
    const Point_map_3& point_map_3, 
    Plane_3& plane) {

    std::vector<std::size_t> indices;
    indices.reserve(item_range.size());
    for (std::size_t i = 0; i < item_range.size(); ++i)
      indices.push_back(i);

    return plane_from_points_3(item_range, point_map_3, indices, plane);
  }

  template<typename Vector_3>
  typename Kernel_traits<Vector_3>::Kernel::FT
  vector_length(const Vector_3& v) {

    using Traits = typename Kernel_traits<Vector_3>::Kernel;
    using FT = typename Traits::FT;
    using Get_sqrt = Get_sqrt<Traits>;
    using Sqrt = typename Get_sqrt::Sqrt;
    const Sqrt sqrt;
    return static_cast<FT>(sqrt(v * v));
  }

  template<
  typename Item_range,
  typename Point_map_2,
  typename Line_2,
  typename Point_2>
  void boundary_points_on_line_2(
    const Item_range& item_range,
    const Point_map_2 point_map_2,
    const std::vector<std::size_t>& indices,
    const Line_2& line,
    Point_2& p,
    Point_2& q) {

    using Traits = typename Kernel_traits<Line_2>::Kernel;
    using FT = typename Traits::FT;
    using Vector_2 = typename Traits::Vector_2;

    FT min_proj_value = max_value<FT>();
    FT max_proj_value = -max_value<FT>();

    const Vector_2 ref_vector = line.to_vector();
    const Point_2& ref_point = get(point_map_2, item_range[indices[0]]);
    
    for (std::size_t i = 0; i < indices.size(); ++i) {
      const Point_2& query = get(point_map_2, item_range[indices[i]]);
      const Point_2 point = line.projection(query);
      
      const Vector_2 curr_vector(ref_point, point);
      const FT value = CGAL::scalar_product(curr_vector, ref_vector);
      
      if (value < min_proj_value) {
        min_proj_value = value;
        p = point; }
      if (value > max_proj_value) {
        max_proj_value = value;
        q = point; }
    }
  }

  template<
  typename Item_range,
  typename Point_map_2,
  typename Line_2>
  typename Kernel_traits<Line_2>::Kernel::FT
  points_squared_length_2(
    const Item_range& item_range,
    const Point_map_2& point_map_2,
    const std::vector<std::size_t>& indices,
    const Line_2& line) {

    using Traits = typename Kernel_traits<Line_2>::Kernel;
    using FT = typename Traits::FT;
    using Point_2 = typename Traits::Point_2;

    Point_2 p, q;
    boundary_points_on_line_2(item_range, point_map_2, indices, line, p, q);
    const FT squared_length = 
    CGAL::squared_distance(line.projection(p), line.projection(q));
    return squared_length;
  }

  template<
  typename Item_range, 
  typename Point_map_2, 
  typename Line_2>
  typename Kernel_traits<Line_2>::Kernel::FT
  line_from_points_2(
    const Item_range& item_range, 
    const Point_map_2& point_map_2, 
    Line_2& line) {

    using Traits = typename Kernel_traits<Line_2>::Kernel;
    using FT = typename Traits::FT;

    using Local_traits 
    = CGAL::Exact_predicates_inexact_constructions_kernel;
		using Local_FT = typename Local_traits::FT;
		using Local_line_2 = typename Local_traits::Line_2;
    using Local_point_2 = typename Local_traits::Point_2;

		CGAL_assertion(item_range.size() > 0);
		std::vector<Local_point_2> points;
    points.reserve(item_range.size());
				
		for (std::size_t i = 0; i < item_range.size(); ++i) {
			const auto& p = get(point_map_2, *(item_range.begin() + i));

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

  template<typename Point>
  typename Kernel_traits<Point>::Kernel::FT
  distance(
    const Point& p, 
    const Point& q) {
      
    using Traits = typename Kernel_traits<Point>::Kernel;
    using FT = typename Traits::FT;
    return static_cast<FT>(
      CGAL::sqrt(CGAL::to_double(CGAL::squared_distance(p, q))));
  }

} // internal
} // Urban_area_processing
} // CGAL

#endif // CGAL_URBAN_AREA_PROCESSING_INTERNAL_UTILS_H
