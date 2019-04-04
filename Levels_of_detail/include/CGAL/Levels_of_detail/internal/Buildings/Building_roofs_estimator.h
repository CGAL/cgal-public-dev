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

#ifndef CGAL_LEVELS_OF_DETAIL_INTERNAL_ROOFS_ESTIMATOR_H
#define CGAL_LEVELS_OF_DETAIL_INTERNAL_ROOFS_ESTIMATOR_H

// STL includes.
#include <vector>

// CGAL includes.
#include <CGAL/assertions.h>

// Internal includes.
#include <CGAL/Levels_of_detail/internal/utils.h>
#include <CGAL/Levels_of_detail/internal/struct.h>

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

  template<
  typename GeomTraits, 
  typename InputRange,
  typename PointMap>
  class Building_roofs_estimator {

  public:
    using Traits = GeomTraits;
    using Input_range = InputRange;
    using Point_map = PointMap;

    using FT = typename Traits::FT;
    using Point_2 = typename Traits::Point_2;
    using Point_3 = typename Traits::Point_3;
    using Plane_3 = typename Traits::Plane_3;
    using Vector_2 = typename Traits::Vector_2;
    using Vector_3 = typename Traits::Vector_3;

    using Indices = std::vector<std::size_t>;
    using Points_3 = std::vector<Point_3>;
    using Approximate_face = internal::Partition_edge_3<Traits>;

    Building_roofs_estimator(
      const Input_range& input_range, 
      const Point_map& point_map,
      const std::vector<Indices>& roof_points_3) :
    m_input_range(input_range),
    m_point_map(point_map),
    m_roof_points_3(roof_points_3)
    { }

    void estimate(
      std::vector<Approximate_face>& roofs) const {

      roofs.clear();      
      if (m_roof_points_3.empty())
        return;

      roofs.reserve(m_roof_points_3.size());
      for (const auto& indices : m_roof_points_3)
        estimate_roof(indices, roofs);
    }

  private:
    const Input_range& m_input_range;
    const Point_map& m_point_map;
    const std::vector<Indices>& m_roof_points_3;

    void estimate_roof(
      const Indices& indices,
      std::vector<Approximate_face>& roofs) const {

      CGAL_assertion(indices.size() > 0);

      Plane_3 plane;
      internal::plane_from_points_3(
        m_input_range, m_point_map, indices, plane);

      Vector_3 m = plane.orthogonal_vector();
      Vector_3 n = Vector_3(FT(0), FT(0), FT(1));
      if (m == -n) m = n;

      FT angle_3d; Vector_3 axis;
      const bool success = internal::compute_angle_and_axis_3(
        m, n, angle_3d, axis);
      if (!success) return;
      const FT angle_deg = angle_3d * FT(180) / static_cast<FT>(CGAL_PI);

      Points_3 points;
      internal::project_on_plane_3(
        m_input_range, m_point_map, indices, plane, points);

			if (angle_deg != FT(0) && angle_deg != FT(180))
				internal::rotate_points_3(angle_3d, axis, points);

      Vector_2 dir;
      internal::estimate_direction_2(points, dir);
      const Vector_2 y_dir = Vector_2(FT(0), FT(1));

      FT angle_2d;
			internal::compute_angle_2(dir, y_dir, angle_2d);

			Point_3 bar;
      internal::compute_barycenter_3(points, bar);
      const Point_2 b = Point_2(bar.x(), bar.y());
      rotate_points(angle_2d, b, points);

			Approximate_face face;
			compute_bounding_box(points, face.polygon);

      rotate_points(-angle_2d, b, face.polygon);
      if (angle_deg != FT(0) && angle_deg != FT(180))
				internal::rotate_points_3(-angle_3d, axis, face.polygon);

			if (!is_valid_boundary(face.polygon)) {
				face.polygon.clear(); return; }
      roofs.push_back(face);
    }

    void rotate_points(
      const FT angle, 
      const Point_2& barycenter, 
      Points_3& points) const {

      Point_2 q; FT z;
      for (auto& p : points) {
        q = Point_2(p.x(), p.y()); z = p.z();
        internal::rotate_point_2(angle, barycenter, q);
        p = Point_3(q.x(), q.y(), z);
      }
    }

		void compute_bounding_box(
      const Points_3& points, Points_3& bbox) const {

			FT minx = internal::max_value<FT>(), miny = internal::max_value<FT>();
			FT maxx = -internal::max_value<FT>(), maxy = -internal::max_value<FT>();

			FT z = FT(0);
			for (const auto& p : points) {
				minx = CGAL::min(minx, p.x()); miny = CGAL::min(miny, p.y());
				maxx = CGAL::max(maxx, p.x()); maxy = CGAL::max(maxy, p.y());
				z += p.z();
			}
			z /= static_cast<FT>(points.size());

      bbox.clear(); bbox.reserve(4);
      bbox.push_back(Point_3(minx, miny, z));
			bbox.push_back(Point_3(maxx, miny, z));
			bbox.push_back(Point_3(maxx, maxy, z));
			bbox.push_back(Point_3(minx, maxy, z));
      CGAL_assertion(bbox.size() == 4);
		}

		bool is_valid_boundary(
      const Points_3& boundary) const {

			if (std::isnan(CGAL::to_double(boundary[0].x())) ||
					std::isnan(CGAL::to_double(boundary[0].y())) ||
					std::isnan(CGAL::to_double(boundary[0].z()))  ) return false;

      if (std::isnan(CGAL::to_double(boundary[1].x())) ||
					std::isnan(CGAL::to_double(boundary[1].y())) ||
					std::isnan(CGAL::to_double(boundary[1].z()))  ) return false;

      if (std::isnan(CGAL::to_double(boundary[2].x())) ||
					std::isnan(CGAL::to_double(boundary[2].y())) ||
					std::isnan(CGAL::to_double(boundary[2].z()))  ) return false;

      if (std::isnan(CGAL::to_double(boundary[3].x())) ||
					std::isnan(CGAL::to_double(boundary[3].y())) ||
					std::isnan(CGAL::to_double(boundary[3].z()))  ) return false;

			if (boundary.size() < 3) 
        return false;
      return true;
		}
  };

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_INTERNAL_ROOFS_ESTIMATOR_H
