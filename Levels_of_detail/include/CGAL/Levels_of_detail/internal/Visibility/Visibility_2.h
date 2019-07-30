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

#ifndef CGAL_LEVELS_OF_DETAIL_INTERNAL_VISIBILITY_2_H
#define CGAL_LEVELS_OF_DETAIL_INTERNAL_VISIBILITY_2_H

// STL includes.
#include <vector>
#include <utility>

// CGAL includes.
#include <CGAL/assertions.h>
#include <CGAL/property_map.h>

// LOD includes.
#include <CGAL/Levels_of_detail/enum.h>

// Internal includes.
#include <CGAL/Levels_of_detail/internal/utils.h>
#include <CGAL/Levels_of_detail/internal/struct.h>

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

template<
typename GeomTraits,
typename PointMap>
class Visibility_2 {

public:
  using Traits = GeomTraits;
  using Point_map = PointMap;

  using FT = typename Traits::FT;
  using Point_2 = typename Traits::Point_2;
  using Triangle_2 = typename Traits::Triangle_2;

  using Points_2 = std::vector<Point_2>;
  using Indices = std::vector<std::size_t>;

  using Triangulation = typename internal::Triangulation<Traits>::Delaunay;
  using Location_type = typename Triangulation::Locate_type;

  using Partition_2 = internal::Partition_2<Traits>;

  Visibility_2(
    const Indices& boundary_points,
    const Indices& interior_points,
    const Point_map& point_map,
    const FT alpha,
    const FT threshold = FT(1) / FT(2)) :
  m_boundary_points(boundary_points),
  m_interior_points(interior_points),
  m_point_map(point_map),
  m_alpha(alpha),
  m_threshold(threshold)
  { }

  void compute(Partition_2& partition_2) const {

    // Create input items.
    std::vector< std::pair<std::size_t, bool> > items;
    items.reserve(m_boundary_points.size() + m_interior_points.size());
    for (const std::size_t idx : m_boundary_points)
      items.push_back(std::make_pair(idx, false));
    for (const std::size_t idx : m_interior_points)
      items.push_back(std::make_pair(idx, false));

    auto pmap = CGAL::Identity_property_map<Point_2>();
    Points_2 points; Indices indices;
    for (auto& face : partition_2.faces) {
      const auto& tri = face.base.delaunay;

      // Find face points.
      points.clear(); indices.clear(); std::size_t idx = 0;
      for (auto& item : items) {
        if (!item.second) {
            
          const auto& p = get(m_point_map, item.first);
          Location_type type; int stub;
          const auto fh = tri.locate(p, type, stub);
          if (type == Triangulation::FACE && !tri.is_infinite(fh)) {
            points.push_back(p);
            indices.push_back(idx);
            item.second = true; ++idx;
          }
        }
      }

      // Compute total face area.
      FT face_area = FT(0);
      for (auto fh = tri.finite_faces_begin(); 
      fh != tri.finite_faces_end(); ++fh) {

        const Triangle_2 triangle = Triangle_2(
          fh->vertex(0)->point(),
          fh->vertex(1)->point(),
          fh->vertex(2)->point());
        face_area += triangle.area();
      }

      // Compute approximate alpha_shapes-based face area.
      FT approx_area = FT(0);
      if (points.size() >= 3)
        approx_area = internal::points_area(points, pmap, indices, m_alpha);

      // Set labels.
      if (face_area > approx_area) face.inside = approx_area / face_area;
      else face.inside = face_area / approx_area;
      face.outside = FT(1) - face.inside;

      if (face.inside > m_threshold)
        face.visibility = Visibility_label::INSIDE;
      else
        face.visibility = Visibility_label::OUTSIDE;
    }
  }

private:
  const Indices& m_boundary_points;
  const Indices& m_interior_points;
  const Point_map& m_point_map;
  const FT m_alpha;
  const FT m_threshold;
};

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_INTERNAL_VISIBILITY_2_H
