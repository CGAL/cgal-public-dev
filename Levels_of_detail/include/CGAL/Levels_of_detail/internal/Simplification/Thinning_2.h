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

#ifndef CGAL_LEVELS_OF_DETAIL_INTERNAL_THINNING_2_H
#define CGAL_LEVELS_OF_DETAIL_INTERNAL_THINNING_2_H

// STL includes.
#include <vector>

// Internal includes.
#include <CGAL/Levels_of_detail/internal/utils.h>

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

  template<
  typename GeomTraits,
  typename NeighborQuery>
  class Thinning_2 {

  public:
    using Traits = GeomTraits;
    using Neighbor_query = NeighborQuery;

    using FT = typename Traits::FT;
    using Point_2 = typename Traits::Point_2;
    using Line_2 = typename Traits::Line_2;

    Thinning_2(
      const Neighbor_query& neighbor_query) :
    m_neighbor_query(neighbor_query)
    { }

    void apply(std::vector<Point_2>& points) const {

      if (points.empty())
        return;

      std::vector<Point_2> tmp;
      tmp.reserve(points.size());
      Line_2 line;
      std::vector<std::size_t> neighbors;

      for (const Point_2& p : points) {
        m_neighbor_query(p, neighbors);
        internal::line_from_points_2(
          neighbors, m_neighbor_query.point_map(), line);
        tmp.push_back(line.projection(p));
      }
      points = tmp;
    }

  private:
    const Neighbor_query& m_neighbor_query;
  };

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_INTERNAL_THINNING_2_H
