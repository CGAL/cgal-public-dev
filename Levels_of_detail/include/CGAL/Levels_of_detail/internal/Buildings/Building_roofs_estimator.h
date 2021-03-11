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
#include <CGAL/property_map.h>

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
    using Points_2 = std::vector<Point_2>;
    using Points_3 = std::vector<Point_3>;
    using Approximate_face = internal::Partition_edge_3<Traits>;
    using Identity_map_2 = CGAL::Identity_property_map<Point_2>;

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

      Points_3 points_3;
      internal::project_on_plane_3(
        m_input_range, m_point_map, indices, plane, points_3);

      Point_3 b3;
      internal::compute_barycenter_3(points_3, b3);

      Points_2 points_2;
      points_2.reserve(points_3.size());
      for (const auto& p3 : points_3)
        points_2.push_back(internal::to_2d(p3, b3, plane));

      Vector_2 dir;
      internal::estimate_direction_2(points_2, dir);
      const Vector_2 y_dir = Vector_2(FT(0), FT(1));

      FT angle_2d;
			internal::compute_angle_2(dir, y_dir, angle_2d);

			Point_2 b2;
      internal::compute_barycenter_2(points_2, b2);
      internal::rotate_points_2(angle_2d, b2, points_2);

			Points_2 bbox;
      Identity_map_2 identity_map_2;
			internal::bounding_box_2(points_2, identity_map_2, bbox);

      Approximate_face face; face.polygon.reserve(4);
      internal::rotate_points_2(-angle_2d, b2, bbox);
      for (const auto& p : bbox)
        face.polygon.push_back(internal::to_3d(p, b3, plane));

			if (!is_valid_boundary(face.polygon)) {
				face.polygon.clear(); return; }
      roofs.push_back(face);
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
