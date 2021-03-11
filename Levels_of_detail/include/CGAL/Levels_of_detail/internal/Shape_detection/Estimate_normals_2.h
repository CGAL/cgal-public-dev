// Copyright (c) 2019 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
//
//
// Author(s)     : Dmitry Anisimov, Simon Giraudot, Pierre Alliez, Florent Lafarge, and Andreas Fabri
//

#ifndef CGAL_LEVELS_OF_DETAIL_INTERNAL_SHAPE_DETECTION_ESTIMATE_NORMALS_2_H
#define CGAL_LEVELS_OF_DETAIL_INTERNAL_SHAPE_DETECTION_ESTIMATE_NORMALS_2_H

// STL includes.
#include <vector>

// CGAL includes.
#include <CGAL/assertions.h>

// Internal includes.
#include <CGAL/Levels_of_detail/internal/utils.h>
#include <CGAL/Levels_of_detail/internal/property_map.h>

// Testing.
#include "../../../../../test/Levels_of_detail/include/Saver.h"

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

  template<
  typename GeomTraits,
  typename InputRange,
  typename PointMap,
  typename NeighborQuery>
  class Estimate_normals_2 {

  public:
    using Traits = GeomTraits;
    using Input_range = InputRange;
    using Point_map = PointMap;
    using Neighbor_query = NeighborQuery;

    using FT = typename Traits::FT;
    using Point_2 = typename Traits::Point_2;
    using Point_3 = typename Traits::Point_3;
    using Vector_2 = typename Traits::Vector_2;
    using Vector_3 = typename Traits::Vector_3;
    using Line_2 = typename Traits::Line_2;

    using Indices = std::vector<std::size_t>;
    using Saver = Saver<Traits>;

    Estimate_normals_2(
      const Input_range& input_range,
      const Neighbor_query& neighbor_query,
      const Point_map& point_map) :
    m_input_range(input_range),
    m_neighbor_query(neighbor_query),
    m_point_map(point_map) {
      CGAL_precondition(m_input_range.size() > 0);
    }

    void get_normals(std::vector<Vector_2>& normals) const {

      normals.clear();
      normals.reserve(m_input_range.size());

      Line_2 line; Vector_2 normal;
      Indices neighbors;
      for (std::size_t i = 0; i < m_input_range.size(); ++i) {

        neighbors.clear();
        m_neighbor_query(i, neighbors);
        internal::line_from_points_2(
          neighbors, m_neighbor_query.point_map(), line);

        normal = line.to_vector();
        normal = normal.perpendicular(CGAL::COUNTERCLOCKWISE);

        const FT normal_length = internal::vector_length(normal);
        CGAL_precondition(normal_length > FT(0));
        normal /= normal_length;
        normals.push_back(normal);
      }
      CGAL_assertion(normals.size() == m_input_range.size());
      save_normals(normals);
    }

  private:
    const Input_range& m_input_range;
    const Neighbor_query& m_neighbor_query;
    const Point_map& m_point_map;

    void save_normals(
      const std::vector<Vector_2>& normals) const {

      Saver saver;
      std::vector<Point_3> points3(m_input_range.size());
      std::vector<Vector_3> normals3(m_input_range.size());
      for (std::size_t i = 0; i < m_input_range.size(); ++i) {
        const Point_2& p = get(m_point_map, *(m_input_range.begin() + i));
        const Vector_2& n = normals[i];

        points3[i]= Point_3(p.x(), p.y(), FT(0));
        normals3[i] = Vector_3(n.x(), n.y(), FT(0));
      }

      saver.export_points_with_normals(
        points3, normals3, "/Users/monet/Documents/gf/lod/logs/buildings/tmp/normals");
    }
  };

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_INTERNAL_SHAPE_DETECTION_ESTIMATE_NORMALS_2_H
