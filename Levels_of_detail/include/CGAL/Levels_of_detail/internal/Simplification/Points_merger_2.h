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

#ifndef CGAL_LEVELS_OF_DETAIL_INTERNAL_POINTS_MERGER_2_H
#define CGAL_LEVELS_OF_DETAIL_INTERNAL_POINTS_MERGER_2_H

// STL includes.
#include <map>
#include <utility>
#include <iostream>
#include <vector>
#include <string>
#include <list>
#include <limits>
#include <set>
#include <algorithm>
#include <iterator>

// CGAL includes.
#include <CGAL/barycenter.h>
#include <CGAL/property_map.h>

// Internal includes.
#include <CGAL/Levels_of_detail/internal/utils.h>
#include <CGAL/Levels_of_detail/internal/struct.h>

// Spatial search.
#include <CGAL/Levels_of_detail/internal/Spatial_search/K_neighbor_query.h>
#include <CGAL/Levels_of_detail/internal/Spatial_search/Sphere_neighbor_query.h>

// Testing.
#include "../../../../../test/Levels_of_detail/include/Saver.h"

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

  template<typename GeomTraits>
  class Points_merger_2 {

  public:
    using Traits = GeomTraits;
    using FT = typename Traits::FT;
    using Point_2 = typename Traits::Point_2;
    using Point_3 = typename Traits::Point_3;
    using Segment_2 = typename Traits::Segment_2;
    using Line_2 = typename Traits::Line_2;

    using Pair = std::pair<Point_2, std::size_t>;
    using Point_map = CGAL::First_of_pair_property_map<Pair>;

    using K_neighbor_query =
      internal::K_neighbor_query<Traits, std::vector<Pair>, Point_map>;
    using Sphere_neighbor_query =
      internal::Sphere_neighbor_query<Traits, std::vector<Pair>, Point_map>;

    using Indices = std::vector<std::size_t>;

    using Triangulation = internal::Triangulation<Traits>;

    Points_merger_2(
      const FT noise_level,
      const FT alpha) :
    m_noise_level(noise_level),
    m_alpha(alpha)
    { }

    template<typename Point_map>
    void merge(
      const Indices& indices,
      const Point_map& point_map,
      const std::vector< std::vector<Segment_2> >& contours,
      Triangulation& triangulation) const {
      
      
    }

  private:
    const FT m_noise_level;
    const FT m_alpha;

    void save_points(
      const std::vector<Point_2>& apoints_2,
      const std::string name) const {

      std::vector<Point_3> points;
      points.reserve(apoints_2.size());
      for (const auto& p : apoints_2)
        points.push_back(Point_3(p.x(), p.y(), FT(0)));

      Saver<Traits> saver;
      saver.export_points(
        points, Color(0, 0, 0), name);
    }
  };

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_INTERNAL_POINTS_MERGER_2_H
