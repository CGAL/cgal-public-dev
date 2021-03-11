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
#include <cmath>

// Boost headers.
#include <boost/mpl/has_xxx.hpp>

// CGAL includes.
#include <CGAL/utils.h>
#include <CGAL/Random.h>
#include <CGAL/assertions.h>
#include <CGAL/property_map.h>
#include <CGAL/number_utils.h>
#include <CGAL/Kernel/global_functions.h>
#include <CGAL/compute_average_spacing.h>
#include <CGAL/linear_least_squares_fitting_2.h>
#include <CGAL/linear_least_squares_fitting_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Barycentric_coordinates_2/Segment_coordinates_2.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Alpha_shape_2.h>
#include <CGAL/Polygon_2_algorithms.h>

namespace CGAL {
namespace Levels_of_detail {
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

  template<typename Point_3>
  bool are_coplanar_3(
    const Point_3& p1, const Point_3& p2,
    const Point_3& p3, const Point_3& p4) {

    using Traits = typename Kernel_traits<Point_3>::Kernel;
    using FT = typename Traits::FT;
    using Vector_3 = typename Traits::Vector_3;

		const Vector_3 v1 = Vector_3(p1, p2);
		const Vector_3 v2 = Vector_3(p1, p3);
		const Vector_3 v3 = Vector_3(p1, p4);

		const Vector_3 v4 = CGAL::cross_product(v2, v3);
		const FT result = CGAL::scalar_product(v1, v4);
		return CGAL::abs(result) < tolerance<FT>();
	}

  template<typename Point_2>
  bool are_equal_points_2(
    const Point_2& p, const Point_2& q) {

    using Traits = typename Kernel_traits<Point_2>::Kernel;
    using FT = typename Traits::FT;

    const FT eps = tolerance<FT>();
    return
    (CGAL::abs(p.x() - q.x()) < eps) &&
    (CGAL::abs(p.y() - q.y()) < eps);
  }

  template<typename Point_3>
  bool are_equal_points_3(
    const Point_3& p, const Point_3& q) {

    using Traits = typename Kernel_traits<Point_3>::Kernel;
    using FT = typename Traits::FT;

    const FT eps = tolerance<FT>();
    return
    (CGAL::abs(p.x() - q.x()) < eps) &&
    (CGAL::abs(p.y() - q.y()) < eps) &&
    (CGAL::abs(p.z() - q.z()) < eps);
  }

  template<typename Point_3>
	bool are_equal_edges_3(
      const Point_3& p1, const Point_3& p2,
      const Point_3& q1, const Point_3& q2) {

    return (
      internal::are_equal_points_3(p1, q1) &&
      internal::are_equal_points_3(p2, q2)) || (
      internal::are_equal_points_3(p1, q2) &&
      internal::are_equal_points_3(p2, q1));
  }

  template<typename Vector>
  void normalize(Vector& v) {

    using Traits = typename Kernel_traits<Vector>::Kernel;
    using FT = typename Traits::FT;

    v /= static_cast<FT>(
      CGAL::sqrt(CGAL::to_double(v.squared_length())));
  }

  template<
  typename Point_3,
  typename Vector_3>
  bool compute_cross_product_3(
    const std::vector<Point_3>& polygon,
    Vector_3& cross) {

    using Traits = typename Kernel_traits<Point_3>::Kernel;
    using FT = typename Traits::FT;

    CGAL_assertion(polygon.size() >= 3);
    const Vector_3 zero = Vector_3(FT(0), FT(0), FT(0));
    for (std::size_t i = 0; i < polygon.size(); ++i) {

      const std::size_t ip = (i + 1) % polygon.size();
      const std::size_t ipp = (i + 2) % polygon.size();

      const Point_3& p1 = polygon[i];
      const Point_3& p2 = polygon[ip];
      const Point_3& p3 = polygon[ipp];

      const Vector_3 v1 = Vector_3(p2, p1);
      const Vector_3 v2 = Vector_3(p2, p3);

      cross = CGAL::cross_product(v1, v2);
      if (!are_equal_points_3(cross, zero))
        return true;
    }
    return false;
  }

  template<
  typename Vector_2,
  typename FT>
  void compute_angle_2(
    const Vector_2& m, const Vector_2& n, FT& angle) {

		const FT cross = CGAL::determinant(m, n);
		const FT dot = CGAL::scalar_product(m, n);
		angle = static_cast<FT>(
      std::atan2(CGAL::to_double(cross), CGAL::to_double(dot)));
	}

  template<typename Vector_2>
  typename Kernel_traits<Vector_2>::Kernel::FT
  unsigned_angle_2(
    const Vector_2& m, const Vector_2& n) {

    using FT = typename Kernel_traits<Vector_2>::Kernel::FT;
		const FT dot = CGAL::scalar_product(m, n);

    const FT lm = static_cast<FT>(
      CGAL::sqrt(CGAL::to_double(m.squared_length())));
    const FT ln = static_cast<FT>(
      CGAL::sqrt(CGAL::to_double(n.squared_length())));

    const FT len = lm * ln;
		return static_cast<FT>(
      std::acos(CGAL::to_double(dot) / CGAL::to_double(len)));
	}

  template<
  typename Point_2,
  typename FT>
  bool is_thin_polygon_2(
    const std::vector<Point_2>& polygon,
    const FT max_angle) {

		using Traits = typename Kernel_traits<Point_2>::Kernel;
    using Vector_2 = typename Traits::Vector_2;

    std::size_t count = 0;
    for (std::size_t i = 0; i < polygon.size(); ++i) {
      const std::size_t ip = (i + 1) % polygon.size();
      const std::size_t ipp = (i + 2) % polygon.size();

      const auto& p0 = polygon[i];
      const auto& p1 = polygon[ip];
      const auto& p2 = polygon[ipp];

      const Vector_2 v1 = Vector_2(p1, p0);
      const Vector_2 v2 = Vector_2(p1, p2);

      FT angle_rad;
      compute_angle_2(v1, v2, angle_rad);
      const FT angle_deg = angle_rad * FT(180) / static_cast<FT>(CGAL_PI);

      if (CGAL::abs(angle_deg) <= max_angle)
        ++count;
    }
    return count >= 2;
	}

  template<
  typename Point_3,
  typename Vector_3>
  bool compute_normal_3(
    const std::vector<Point_3>& polygon,
    Vector_3& normal) {

    CGAL_assertion(polygon.size() >= 3);
    if (polygon.size() < 3)
      return false;

    const bool success =
      compute_cross_product_3(polygon, normal);
    if (success) {
      normalize(normal); return true;
    } return false;
  }

  template<
  typename Vector_3,
  typename FT>
  bool compute_angle_and_axis_3(
    const Vector_3& m, const Vector_3& n,
    FT& angle, Vector_3& axis) {

		const auto cross = CGAL::cross_product(m, n);
		const FT length = static_cast<FT>(
      CGAL::sqrt(CGAL::to_double(cross.squared_length())));
		const FT dot = CGAL::scalar_product(m, n);

		angle = static_cast<FT>(std::atan2(
      CGAL::to_double(length), CGAL::to_double(dot)));

		const FT angle_deg =
      angle * FT(180) / static_cast<FT>(CGAL_PI);
		if (angle_deg == FT(0) || angle_deg == FT(180))
			return true;

		if (length == FT(0)) {
      std::cerr << "Error: length = 0" << std::endl;
      exit(EXIT_FAILURE); }

		CGAL_assertion(length > FT(0));
		axis = cross / length;

    const FT half_pi = static_cast<FT>(CGAL_PI) / FT(2);
    if (angle > half_pi) {
      angle = static_cast<FT>(CGAL_PI) - angle;
      axis = -axis;
    }
		return true;
	}

  template<
  typename Vector_3,
  typename FT>
  void compute_angle_3_deg(
    const Vector_3& m, const Vector_3& n,
    FT& angle_deg) {

		const auto cross = CGAL::cross_product(m, n);
		const FT length = static_cast<FT>(
      CGAL::sqrt(CGAL::to_double(cross.squared_length())));
		const FT dot = CGAL::scalar_product(m, n);

		FT angle_rad = static_cast<FT>(std::atan2(
      CGAL::to_double(length), CGAL::to_double(dot)));

    const FT half_pi = static_cast<FT>(CGAL_PI) / FT(2);
    if (angle_rad > half_pi)
        angle_rad = static_cast<FT>(CGAL_PI) - angle_rad;
		angle_deg = angle_rad * FT(180) / static_cast<FT>(CGAL_PI);
  }

  template<typename Point_2>
  void compute_barycenter_2(
    const std::vector<Point_2>& points,
    Point_2& b) {

    using Traits = typename Kernel_traits<Point_2>::Kernel;
    using FT = typename Traits::FT;

    CGAL_assertion(points.size() > 0);
    FT x = FT(0), y = FT(0);
    for (const auto& p : points) {
      x += p.x();
      y += p.y();
    }
    x /= static_cast<FT>(points.size());
    y /= static_cast<FT>(points.size());
    b = Point_2(x, y);
  }

  template<typename Point_3>
  void compute_barycenter_3(
    const std::vector<Point_3>& points,
    Point_3& b) {

    using Traits = typename Kernel_traits<Point_3>::Kernel;
    using FT = typename Traits::FT;

    CGAL_assertion(points.size() > 0);
    FT x = FT(0), y = FT(0), z = FT(0);
    for (const auto& p : points) {
      x += p.x();
      y += p.y();
      z += p.z();
    }
    x /= static_cast<FT>(points.size());
    y /= static_cast<FT>(points.size());
    z /= static_cast<FT>(points.size());
    b = Point_3(x, y, z);
  }

  template<typename Point_3>
  void compute_barycenter_3(
    const std::vector< std::vector<Point_3> >& polygons,
    const std::vector<std::size_t>& indices,
    Point_3& b) {

    using Traits = typename Kernel_traits<Point_3>::Kernel;
    using FT = typename Traits::FT;

    CGAL_assertion(polygons.size() > 0 && indices.size() > 0);
    Point_3 sub_b; FT x = FT(0), y = FT(0), z = FT(0);
    for (const std::size_t idx : indices) {
      CGAL_assertion(idx >= 0 && idx < polygons.size());
      const auto& polygon = polygons[idx];
      compute_barycenter_3(polygon, sub_b);

      x += sub_b.x();
      y += sub_b.y();
      z += sub_b.z();
    }
    x /= static_cast<FT>(indices.size());
    y /= static_cast<FT>(indices.size());
    z /= static_cast<FT>(indices.size());
    b = Point_3(x, y, z);
  }

  template<
  typename FT,
  typename Point_2>
  void scale_polygon_2(
    const FT scale,
    std::vector<Point_2>& polygon) {

    Point_2 b;
    compute_barycenter_2(polygon, b);

    for (auto& p : polygon) {
      const FT x = (p.x() - b.x()) * scale + b.x();
      const FT y = (p.y() - b.y()) * scale + b.y();
      p = Point_2(x, y);
    }
  }

  template<
  typename FT,
  typename Point_3>
  void scale_polygon_3(
    const FT scale,
    const FT z_extender,
    std::vector<Point_3>& polygon) {

    Point_3 b;
    compute_barycenter_3(polygon, b);

    for (auto& p : polygon) {
      const FT x = (p.x() - b.x()) * scale + b.x();
      const FT y = (p.y() - b.y()) * scale + b.y();
      const FT z = (p.z() - b.z()) * scale * z_extender + b.z();
      p = Point_3(x, y, z);
    }
  }

  template<
  typename FT,
  typename Vector_3,
  typename Point_3>
  void rotate_point_3(
    const FT angle,
    const Vector_3& axis,
    Point_3& p) {

		const double tmp_angle = CGAL::to_double(angle);
		const FT c = static_cast<FT>(std::cos(tmp_angle));
		const FT s = static_cast<FT>(std::sin(tmp_angle));
		const FT C = FT(1) - c;

		const FT x = axis.x();
		const FT y = axis.y();
		const FT z = axis.z();

		p = Point_3(
      (x * x * C + c)     * p.x() + (x * y * C - z * s) * p.y() + (x * z * C + y * s) * p.z(),
			(y * x * C + z * s) * p.x() + (y * y * C + c)     * p.y() + (y * z * C - x * s) * p.z(),
			(z * x * C - y * s) * p.x() + (z * y * C + x * s) * p.y() + (z * z * C + c)     * p.z());
	}

  template<
  typename FT,
  typename Vector_3,
  typename Point_3>
	void rotate_points_3(
    const FT angle,
    const Vector_3& axis,
    std::vector<Point_3>& points) {

		for (auto& p : points)
			rotate_point_3(angle, axis, p);
	}

  template<
  typename FT,
  typename Point_2>
  void rotate_point_2(
    const FT angle,
    const Point_2& barycenter,
    Point_2& p) {

		FT x = p.x(); x -= barycenter.x();
		FT y = p.y(); y -= barycenter.y();

    p = Point_2(x, y);
    const double tmp_angle = CGAL::to_double(angle);
    const FT c = static_cast<FT>(std::cos(tmp_angle));
		const FT s = static_cast<FT>(std::sin(tmp_angle));

		x = p.x() * c - p.y() * s; x += barycenter.x();
		y = p.y() * c + p.x() * s; y += barycenter.y();
		p = Point_2(x, y);
	}

  template<
  typename FT,
  typename Point_2>
	void rotate_points_2(
    const FT angle,
    const Point_2& barycenter,
    std::vector<Point_2>& points) {

		for (auto& p : points)
			rotate_point_2(angle, barycenter, p);
	}

  template<
  typename Point_3,
  typename FT,
  typename Vector_3>
  void rotate_polygon_3(
    const std::vector<Point_3>& polygon,
    const FT angle,
    const Vector_3& axis,
    const Point_3& b,
    std::vector<Point_3>& rotated) {

    if (angle == FT(0)) {
      rotated = polygon; return;
    }

    rotated.clear();
    rotated.reserve(polygon.size());
    Point_3 q;
    for (const auto& p : polygon) {
      q = Point_3(p.x() - b.x(), p.y() - b.y(), p.z() - b.z());
      rotate_point_3(angle, axis, q);
      rotated.push_back(Point_3(q.x() + b.x(), q.y() + b.y(), q.z() + b.z()));
    }
    CGAL_assertion(rotated.size() == polygon.size());
  }

  template<
  typename FT,
  typename Vector_3,
  typename Point_3>
  void rotate_polygon_3(
    const FT angle,
    const Vector_3& axis,
    const Point_3& b,
    std::vector<Point_3>& polygon) {

    if (angle == FT(0)) return;
    Point_3 q;
    for (auto& p : polygon) {
      q = Point_3(p.x() - b.x(), p.y() - b.y(), p.z() - b.z());
      rotate_point_3(angle, axis, q);
      p = Point_3(q.x() + b.x(), q.y() + b.y(), q.z() + b.z());
    }
  }

  template<
  typename Point_3,
  typename FT,
  typename Vector_3>
  void rotate_polygons_3(
    const std::vector< std::vector<Point_3> >& polygons,
    const std::vector<std::size_t>& indices,
    const FT angle,
    const Vector_3& axis,
    const Point_3& b,
    std::vector< std::vector<Point_3> >& rotated) {

    rotated.clear();
    rotated.resize(indices.size());
    for (std::size_t i = 0; i < indices.size(); ++i)
      rotate_polygon_3(polygons[indices[i]], angle, axis, b, rotated[i]);
  }

  template<
  typename Point_3,
  typename FT>
  bool is_vertical_polygon(
    const std::vector<Point_3>& polygon,
    const FT angle_threshold) {

    using Traits = typename Kernel_traits<Point_3>::Kernel;
    using Vector_3 = typename Traits::Vector_3;

		Vector_3 m;
		const bool success = compute_normal_3(polygon, m);
    if (!success) return true;
		const Vector_3 n = Vector_3(FT(0), FT(0), FT(1));
    if (m == -n) m = n;

    FT angle_deg;
    compute_angle_3_deg(m, n, angle_deg);
    const FT angle_diff = CGAL::abs(FT(90) - CGAL::abs(angle_deg));
    if (angle_diff < angle_threshold)
      return true;
    return false;
  }

  template<
  typename Line_3,
  typename Plane_3>
  typename Kernel_traits<Line_3>::Kernel::FT
  intersect_3(
    const Line_3& line, const Plane_3& plane) {

    using Traits = typename Kernel_traits<Line_3>::Kernel;
    using FT = typename Traits::FT;
    using Point_3 = typename Traits::Point_3;
    using Intersect_3 = typename Traits::Intersect_3;

		typename CGAL::cpp11::result_of<Intersect_3(Line_3, Plane_3)>::type result
    = CGAL::intersection(line, plane);
    if (result) {
      if (const Line_3* tmp = boost::get<Line_3>(&*result))
        return max_value<FT>();
      else {
        const Point_3* point = boost::get<Point_3>(&*result);
				return (*point).z();
      }
    }
    return max_value<FT>();
  }

  template<
  typename Point_3,
  typename Point_2>
  void polygon_3_to_polygon_2(
    const std::vector<Point_3>& poly_3,
    std::vector<Point_2>& poly_2) {

    poly_2.clear();
    poly_2.reserve(poly_3.size());
    for (const auto& p : poly_3)
      poly_2.push_back(Point_2(p.x(), p.y()));
    CGAL_assertion(poly_2.size() == poly_3.size());
  }

  template<
  typename Face_handle,
  typename Traits>
  typename Traits::FT min_distance_from_edge_2(
    const Face_handle& fh) {

    using FT = typename Traits::FT;
    using Point_2 = typename Traits::Point_2;
    using Triangle_2 = typename Traits::Triangle_2;

    const Point_2 b = CGAL::barycenter(
      fh->vertex(0)->point(), FT(1),
      fh->vertex(1)->point(), FT(1),
      fh->vertex(2)->point(), FT(1));

    const Triangle_2 triangle = Triangle_2(
      fh->vertex(0)->point(),
      fh->vertex(1)->point(),
      fh->vertex(2)->point());

    FT min_dist = max_value<FT>();
    for (std::size_t i = 0; i < 3; ++i) {
      const std::size_t ip = (i + 1) % 3;
      const FT x = (triangle[i].x() + triangle[ip].x()) / FT(2);
      const FT y = (triangle[i].y() + triangle[ip].y()) / FT(2);
      const Point_2 p = Point_2(x, y);
      const FT dist = distance(p, b);
      min_dist = CGAL::min(min_dist, dist);
    }
    return min_dist;
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
  typename Point_2,
  typename Triangle_2,
  typename FT>
  bool is_near_boundary_2(
    const Point_2& query,
    const Triangle_2& triangle,
    const FT tol) {

    using Traits = typename Kernel_traits<Point_2>::Kernel;
    using Vector_2 = typename Traits::Vector_2;

    std::vector<Vector_2> s(3);
    std::vector<FT> r(3);
    for (std::size_t i = 0; i < 3; ++i) {
      s[i] = triangle[i] - query;
      r[i] = s[i].squared_length();
      if (CGAL::abs(r[i]) < tol * tol)
        return true;
    }

    for (std::size_t i = 0; i < 3; ++i) {
      const std::size_t ip = (i + 1) % 3;
      const FT A = CGAL::determinant(s[i], s[ip]) / FT(2);
      const FT D = CGAL::scalar_product(s[i], s[ip]);
      if (CGAL::abs(A) < tol && D < FT(0))
        return true;
    }
    return false;
  }

  template<
  typename Point_2,
  typename Triangle_2,
  typename FT>
  bool is_within_triangle_2(
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
      const auto res = Barycentric_coordinates::
      compute_segment_coordinates_2(p1, p2, projected, traits);

      const FT squared_tolerance = tol * tol;
      const FT bval = -FT(1) / FT(5);
      const FT tval = FT(6) / FT(5);

      if (res[0] > bval && res[1] > bval &&
          res[0] < tval && res[1] < tval &&
          squared_distance < squared_tolerance)
        return true;
    }
    return false;
  }

  template<typename Point_2>
  bool is_inside_polygon_2(
    const Point_2& query,
    const std::vector<Point_2>& polygon) {

    using Traits = typename Kernel_traits<Point_2>::Kernel;
    using FT = typename Traits::FT;
    using Triangle_2 = typename Traits::Triangle_2;

    CGAL_assertion(polygon.size() > 0);
    if (!CGAL::is_simple_2(polygon.begin(), polygon.end()))
      return false;

    const Point_2& ref = polygon[0];
    for (std::size_t i = 1; i < polygon.size() - 1; ++i) {
      const std::size_t ip = i + 1;

      const Point_2& p1 = ref;
      const Point_2& p2 = polygon[i];
      const Point_2& p3 = polygon[ip];
      const Triangle_2 triangle = Triangle_2(p1, p2, p3);

      if (is_within_triangle_2(query, triangle, tolerance<FT>()))
        return true;
    }
    return false;
  }

  template<
  typename Point_3,
  typename FT>
  FT intersect_with_polygon_3(
    const Point_3& p,
    const std::vector<Point_3>& polygon,
    const FT min_z, const FT max_z) {

    using Traits = typename Kernel_traits<Point_3>::Kernel;
    using Vector_3 = typename Traits::Vector_3;
    using Line_3 = typename Traits::Line_3;
    using Plane_3 = typename Traits::Plane_3;

    const Point_3 p1 = Point_3(p.x(), p.y(), min_z - FT(10));
    const Point_3 p2 = Point_3(p.x(), p.y(), max_z + FT(10));
    const Line_3 line = Line_3(p1, p2);

    Vector_3 m;
    const bool success = compute_normal_3(polygon, m);
    CGAL_assertion(success);
    if (!success) return max_value<FT>();
    Point_3 b;
    compute_barycenter_3(polygon, b);
    const Plane_3 plane = Plane_3(b, m);
    return intersect_3(line, plane);
  }

  std::size_t size_t_rand(
    Random& random,
    const std::size_t maxv)  {

    return static_cast<std::size_t>(random.get_int(0, maxv));
  }

  template<
  typename FT,
  typename Point_3>
  void perturb_point_inside_disc(
    const FT disc_radius,
    const std::size_t num_points_in_disc,
    Random& random,
    Point_3& p) {

    using Local_traits = CGAL::Exact_predicates_inexact_constructions_kernel;
		using Local_FT = Local_traits::FT;
    using Local_point_2 = Local_traits::Point_2;
    using Point_creator_2 = Creator_uniform_2<Local_FT, Local_point_2>;

    std::vector<Local_point_2> points_2;
    points_2.reserve(num_points_in_disc);

    CGAL::Random_points_in_disc_2<Local_point_2, Point_creator_2>
    random_points_2(CGAL::to_double(disc_radius));

    CGAL::cpp11::copy_n(random_points_2, num_points_in_disc,
      std::back_inserter(points_2));

    const std::size_t rand_index = size_t_rand(random, num_points_in_disc - 1);
    const Local_point_2& q = points_2[rand_index];

    const FT qx = static_cast<FT>(q.x());
    const FT qy = static_cast<FT>(q.y());
    p = Point_3(p.x() + qx, p.y() + qy, p.z());
  }

  template<
  typename FT,
  typename Point_3>
  void perturb_polygon_vertices_3(
    const FT disc_radius,
    const std::size_t num_points_in_disc,
    Random& random,
    std::vector<Point_3>& polygon) {

    using Traits = typename Kernel_traits<Point_3>::Kernel;
    using Vector_3 = typename Traits::Vector_3;

    Vector_3 m;
    bool success = compute_normal_3(polygon, m);
    if (!success) return;
    const Vector_3 n = Vector_3(FT(0), FT(0), FT(1));

    FT angle_3d; Vector_3 axis;
    success = compute_angle_and_axis_3(m, n, angle_3d, axis);
    if (!success) return;

    Point_3 b;
    compute_barycenter_3(polygon, b);

    const FT angle_deg = angle_3d * FT(180) / static_cast<FT>(CGAL_PI);
    if (angle_deg != FT(0) && angle_deg != FT(180))
      rotate_polygon_3(angle_3d, axis, b, polygon);

    for (std::size_t i = 0; i < polygon.size(); ++i)
      perturb_point_inside_disc(
        disc_radius, num_points_in_disc, random, polygon[i]);

    if (angle_deg != FT(0) && angle_deg != FT(180))
      rotate_polygon_3(-angle_3d, axis, b, polygon);
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

  template<typename Point_2>
  Point_2 middle_point_2(const Point_2& p, const Point_2& q) {
    using Traits = typename Kernel_traits<Point_2>::Kernel;
    using FT = typename Traits::FT;
    return Point_2((p.x() + q.x()) / FT(2), (p.y() + q.y()) / FT(2));
  }

  template<typename Segment_2>
  typename Kernel_traits<Segment_2>::Kernel::Vector_2
  compute_direction(const Segment_2& segment) {
    using Traits = typename Kernel_traits<Segment_2>::Kernel;
    using FT = typename Traits::FT;
    using Vector = typename Traits::Vector_2;

    Vector v = segment.to_vector();
    if (v.y() < FT(0) || (v.y() == FT(0) && v.x() < FT(0)))
      v = -v;
    normalize(v);
    return v;
  }

  template<typename Vector>
  typename Kernel_traits<Vector>::Kernel::FT
  compute_orientation(const Vector& v) {
    using Traits = typename Kernel_traits<Vector>::Kernel;
    using FT = typename Traits::FT;

    const FT atan = static_cast<FT>(std::atan2(CGAL::to_double(v.y()), CGAL::to_double(v.x())));
    FT orientation = atan * FT(180) / static_cast<FT>(CGAL_PI);
    if (orientation < FT(0))
      orientation += FT(180);
    return orientation;
  }

  template<typename Point_2, typename FT>
  Point_2 transform_coordinates(const Point_2 & barycentre, const Point_2 & frame_origin, const FT angle) {

    const FT cos_val = static_cast<FT>(cos(CGAL_PI * CGAL::to_double(angle) / 180.0));
    const FT sin_val = static_cast<FT>(sin(CGAL_PI * CGAL::to_double(angle) / 180.0));

    const FT x = (barycentre.x() - frame_origin.x()) * cos_val + (barycentre.y() - frame_origin.y()) * sin_val;
    const FT y = (barycentre.y() - frame_origin.y()) * cos_val - (barycentre.x() - frame_origin.x()) * sin_val;

    return Point_2(x, y);
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

  template<
  typename Point_3,
  typename Plane_3>
  typename Kernel_traits<Point_3>::Kernel::Point_2
  to_2d(
    const Point_3& p,
    const Point_3& centroid,
    const Plane_3& plane) {

    using Traits = typename Kernel_traits<Point_3>::Kernel;
    using FT = typename Traits::FT;
    using Point_2 = typename Traits::Point_2;
    using Vector_3 = typename Traits::Vector_3;

    const Vector_3 base1 = plane.base1() / static_cast<FT>(CGAL::sqrt(
      CGAL::to_double(plane.base1() * plane.base1())));
    const Vector_3 base2 = plane.base2() / static_cast<FT>(CGAL::sqrt(
      CGAL::to_double(plane.base2() * plane.base2())));

    const Vector_3 v(centroid, p);
    return Point_2(v * base1, v * base2);
  }

  template<
  typename Point_2,
  typename Point_3,
  typename Plane_3>
  Point_3 to_3d(
    const Point_2& p,
    const Point_3& centroid,
    const Plane_3& plane) {

    using Traits = typename Kernel_traits<Point_3>::Kernel;
    using FT = typename Traits::FT;
    using Vector_3 = typename Traits::Vector_3;

    const Vector_3 base1 = plane.base1() / static_cast<FT>(CGAL::sqrt(
      CGAL::to_double(plane.base1() * plane.base1())));
    const Vector_3 base2 = plane.base2() / static_cast<FT>(CGAL::sqrt(
      CGAL::to_double(plane.base2() * plane.base2())));

    return centroid + p.x() * base1 + p.y() * base2;
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
    Point_2 &p,
    Point_2 &q) {

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
  typename Point_map,
  typename FT>
  FT points_area(
    const Item_range& item_range,
    const Point_map& point_map,
    const std::vector<std::size_t>& indices,
    const FT alpha) {

    using Traits = typename Kernel_traits<
    typename boost::property_traits<Point_map>::value_type>::Kernel;
    using Point_2 = typename Traits::Point_2;
    using Triangle_2 = typename Traits::Triangle_2;

    using Vb = CGAL::Alpha_shape_vertex_base_2<Traits>;
    using Fb = CGAL::Alpha_shape_face_base_2<Traits>;
    using Tds = CGAL::Triangulation_data_structure_2<Vb, Fb>;
    using Triangulation_2 = CGAL::Delaunay_triangulation_2<Traits, Tds>;
    using Alpha_shape_2 = CGAL::Alpha_shape_2<Triangulation_2>;

    Triangulation_2 triangulation;
    for (const std::size_t idx : indices) {
      const auto& p = get(point_map, *(item_range.begin() + idx));
      triangulation.insert(Point_2(p.x(), p.y()));
    }

    FT total_area = FT(0);
    Alpha_shape_2 alpha_shape(triangulation, alpha, Alpha_shape_2::GENERAL);

    for (auto fh = alpha_shape.finite_faces_begin();
      fh != alpha_shape.finite_faces_end(); ++fh) {

      const auto type = alpha_shape.classify(fh);
      if (type == Alpha_shape_2::INTERIOR) {

        const auto& p1 = fh->vertex(0)->point();
        const auto& p2 = fh->vertex(1)->point();
        const auto& p3 = fh->vertex(2)->point();

        const Triangle_2 triangle = Triangle_2(p1, p2, p3);
        total_area += triangle.area();
      }
    }
    return total_area;
  }

  template<
  typename Triangle_2,
  typename Point_2>
  void random_point_in_triangle_2(
    const Triangle_2& triangle,
    Point_2& point) {

    using Traits = typename Kernel_traits<Triangle_2>::Kernel;
    using FT = typename Traits::FT;
    using Vector_2 = typename Traits::Vector_2;

    const Vector_2 v(triangle[0], triangle[1]);
    const Vector_2 w(triangle[0], triangle[2]);
    point = triangle[0];
    double rv = rand() / static_cast<double>(RAND_MAX);
    double rw = rand() / static_cast<double>(RAND_MAX);

    if (rv + rw > 1.0) {
      rv = 1.0 - rv;
      rw = 1.0 - rw;
    }

    const FT bv = static_cast<FT>(rv);
    const FT bw = static_cast<FT>(rw);
    point += bv * v;
    point += bw * w;
  }

  template<
  typename Item_range,
  typename Point_map_3,
  typename Line_3>
  typename Kernel_traits<Line_3>::Kernel::FT
  line_from_points_3(
    const Item_range& item_range,
    const Point_map_3& point_map_3,
    const std::vector<std::size_t>& indices,
    Line_3& line) {

    using Traits = typename Kernel_traits<Line_3>::Kernel;
    using FT = typename Traits::FT;
    using Point_3 = typename Traits::Point_3;
    using Vector_3 = typename Traits::Vector_3;

    using Local_traits
    = CGAL::Exact_predicates_inexact_constructions_kernel;
		using Local_FT = typename Local_traits::FT;
		using Local_line_3 = typename Local_traits::Line_3;
    using Local_point_3 = typename Local_traits::Point_3;

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

		Local_line_3 fitted_line;
    Local_point_3 fitted_centroid;

		const FT quality = static_cast<FT>(
      CGAL::linear_least_squares_fitting_3(
        points.begin(), points.end(),
        fitted_line, fitted_centroid,
        CGAL::Dimension_tag<0>()));

    const auto a = fitted_line.point(0);
    const auto b = fitted_line.to_vector();

    const Point_3 p = Point_3(
      static_cast<FT>(a.x()),
      static_cast<FT>(a.y()),
      static_cast<FT>(a.z()));
    const Vector_3 v = Vector_3(
      static_cast<FT>(b.x()),
      static_cast<FT>(b.y()),
      static_cast<FT>(b.z()));

		line = Line_3(p, v);
    return quality;
  }

  template<
  typename Item_range,
  typename Point_map_3,
  typename Plane_3,
  typename Point_3>
  void project_on_plane_3(
    const Item_range& item_range,
    const Point_map_3& point_map_3,
    const std::vector<std::size_t>& indices,
    const Plane_3& plane,
    std::vector<Point_3>& points) {

    CGAL_assertion(indices.size() > 0);
    points.clear();
    points.reserve(indices.size());

    for (const std::size_t idx : indices) {
      const auto& p = get(point_map_3, *(item_range.begin() + idx));
			points.push_back(plane.projection(p));
    }
    CGAL_assertion(points.size() == indices.size());
  }

  template<
  typename Item_range,
  typename Point_map_2,
  typename Line_2,
  typename Point_2>
  void project_on_line_2(
    const Item_range& item_range,
    const Point_map_2& point_map_2,
    const std::vector<std::size_t>& indices,
    const Line_2& line,
    std::vector<Point_2>& points) {

    CGAL_assertion(indices.size() > 0);
    points.clear();
    points.reserve(indices.size());

    for (const std::size_t idx : indices) {
      const auto& p = get(point_map_2, *(item_range.begin() + idx));
			points.push_back(line.projection(p));
    }
    CGAL_assertion(points.size() == indices.size());
  }

  template<
  typename Item_range,
  typename Point_map_3>
  typename Kernel_traits<
  typename boost::property_traits<Point_map_3>::value_type>::Kernel::FT
  average_distance_to_line_3(
    const Item_range& item_range,
    const Point_map_3& point_map_3,
    const std::vector<std::size_t>& indices) {

    using Traits = typename Kernel_traits<
    typename boost::property_traits<Point_map_3>::value_type>::Kernel;
    using FT = typename Traits::FT;
    using Line_3 = typename Traits::Line_3;

    Line_3 line;
    line_from_points_3(item_range, point_map_3, indices, line);

    FT dist = FT(0);
    for (const std::size_t idx : indices) {
      const auto& p = get(point_map_3, *(item_range.begin() + idx));
      const auto& q = line.projection(p);
      dist += distance(p, q);
    }
    CGAL_assertion(indices.size() > 0);
    dist /= static_cast<FT>(indices.size());
    return dist;
  }

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

    using Local_traits = Exact_predicates_inexact_constructions_kernel;
    using Local_point_3 = typename Local_traits::Point_3;

  public:
    std::size_t operator()(const Point& point) {
      const Local_point_3 p = Local_point_3(
        CGAL::to_double(point.x()),
        CGAL::to_double(point.y()),
        CGAL::to_double(point.z()));

      const auto pair = m_indices.insert(
        std::make_pair(
          p,
          m_indices.size()));
      const auto& item = pair.first;
      const std::size_t idx = item->second;
      return idx;
    }
    void clear() { m_indices.clear(); }

  private:
    std::map<Local_point_3, std::size_t> m_indices;
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
    typename CGAL::cpp11::result_of<Intersect_3(Line_3, Plane_3)>::type result
    = CGAL::intersection(line, plane);
    if (result)
      if (const Point_3* p = boost::get<Point_3>(&*result))
        return *p;

    /*
    std::cerr <<
      "Error (position_on_plane): cannot compute the 3D position!"
    << std::endl; */
    return Point_3(FT(0), FT(0), FT(0));
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

  template<
  typename Point_2,
  typename Vector_2>
  void estimate_direction_2(
    const std::vector<Point_2>& points,
    Vector_2& direction) {

    using Traits = typename Kernel_traits<Point_2>::Kernel;
    using Line_2 = typename Traits::Line_2;

    Line_2 line;
    CGAL::Identity_property_map<Point_2> identity_map;
    line_from_points_2(points, identity_map, line);
    direction = line.to_vector();
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
      const auto& p = get(point_map_2, *(item_range.begin() + i));
      minx = CGAL::min(minx, p.x()); miny = CGAL::min(miny, p.y());
      maxx = CGAL::max(maxx, p.x()); maxy = CGAL::max(maxy, p.y());
    }

    bbox.clear(); bbox.reserve(4);
    bbox.push_back(Point_2(minx, miny)); bbox.push_back(Point_2(maxx, miny));
    bbox.push_back(Point_2(maxx, maxy)); bbox.push_back(Point_2(minx, maxy));
    CGAL_assertion(bbox.size() == 4);
  }

  template<typename Point_2>
  void translate_point_2(const Point_2& tr, Point_2& p) {
    using Traits = typename Kernel_traits<Point_2>::Kernel;
    using FT = typename Traits::FT;

    const FT x = p.x() - tr.x();
    const FT y = p.y() - tr.y();
    p = Point_2(x, y);
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

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_INTERNAL_UTILS_H
