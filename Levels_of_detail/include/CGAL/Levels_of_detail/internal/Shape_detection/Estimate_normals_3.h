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

#ifndef CGAL_LEVELS_OF_DETAIL_INTERNAL_SHAPE_DETECTION_ESTIMATE_NORMALS_3_H
#define CGAL_LEVELS_OF_DETAIL_INTERNAL_SHAPE_DETECTION_ESTIMATE_NORMALS_3_H

// STL includes.
#include <vector>

// CGAL includes.
#include <CGAL/assertions.h>
#include <CGAL/number_utils.h>

// Internal includes.
#include <CGAL/Levels_of_detail/internal/utils.h>

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

  template<
  typename GeomTraits,
  typename InputRange,
  typename PointMap, 
  typename NeighborQuery>
  class Estimate_normals_3 {

  public:
    using Traits = GeomTraits;
    using Input_range = InputRange;
    using Point_map = PointMap;
    using Neighbor_query = NeighborQuery;

    using FT = typename Traits::FT;
    using Vector_3 = typename Traits::Vector_3;
    using Plane_3 = typename Traits::Plane_3;

    Estimate_normals_3(
      const Input_range& input_range,
      const Neighbor_query& neighbor_query,
      const Point_map point_map = Point_map()) : 
    m_input_range(input_range),
    m_neighbor_query(neighbor_query),
    m_point_map(point_map) {
      
      CGAL_precondition(m_input_range.size() > 0);
      estimate_normals();
    }

    const std::vector<Vector_3>& normals() const {
      CGAL_precondition(m_normals.size() == m_input_range.size());
      return m_normals;
    }

  private:
    const Input_range& m_input_range;
    const Neighbor_query& m_neighbor_query;
    const Point_map m_point_map;

    std::vector<Vector_3> m_normals;

    void estimate_normals() {
  
      m_normals.clear();
      m_normals.resize(m_input_range.size());
                  
      Plane_3 plane;
      std::vector<std::size_t> neighbors;
      for (std::size_t i = 0; i < m_input_range.size(); ++i) {

        m_neighbor_query(i, neighbors);
        internal::plane_from_points_3(
          m_input_range, m_point_map, neighbors, plane);
        m_normals[i] = plane.orthogonal_vector();

        const FT normal_length = static_cast<FT>(
          CGAL::sqrt(CGAL::to_double(m_normals[i] * m_normals[i])));

        CGAL_precondition(normal_length > FT(0));
        m_normals[i] /= normal_length;
      }
    }
  };

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_INTERNAL_SHAPE_DETECTION_ESTIMATE_NORMALS_3_H
