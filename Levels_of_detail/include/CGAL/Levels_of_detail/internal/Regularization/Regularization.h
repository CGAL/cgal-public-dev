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

#ifndef CGAL_LEVELS_OF_DETAIL_INTERNAL_REGULARIZATION_H
#define CGAL_LEVELS_OF_DETAIL_INTERNAL_REGULARIZATION_H

#include <CGAL/license/Levels_of_detail.h>

#include <CGAL/Levels_of_detail/internal/Regularization/Shape_regularization.h>
#include <CGAL/Levels_of_detail/internal/Regularization/Angle_regularization_2.h>
#include <CGAL/Levels_of_detail/internal/Regularization/Ordinate_regularization_2.h>
#include <CGAL/Levels_of_detail/internal/Regularization/Delaunay_neighbor_query_2.h>

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

  template<typename GeomTraits>
  class Regularization {

  public:
    using Traits = GeomTraits;

    using FT = typename Traits::FT;
    using Segment_2 = typename Traits::Segment_2;

    using Indices = std::vector<std::size_t>;
    using Segment_range = std::vector<Segment_2>;
    using Segment_map = CGAL::Identity_property_map<Segment_2>;

    using Neighbor_query = 
    internal::Delaunay_neighbor_query_2<Traits, Segment_range, Segment_map>;
    using RT_angles = 
    internal::Angle_regularization_2<Traits, Segment_range, Segment_map>;
    using RT_ordinates = 
    internal::Ordinate_regularization_2<Traits, Segment_range, Segment_map>;

    using SR_angles = internal::Shape_regularization
    <Traits, Segment_range, Neighbor_query, RT_angles>;
    using SR_ordinates = internal::Shape_regularization
    <Traits, Segment_range, Neighbor_query, RT_ordinates>;

    Regularization(
      Segment_range& segments):
    m_segments(segments) 
    { }

    void regularize_angles(const FT angle_bound) {
      std::cout << "Angle bound: " << angle_bound << std::endl;
      
      Neighbor_query neighbor_query(m_segments);

      std::vector<std::size_t> group;
      group.resize(m_segments.size());
      std::iota(group.begin(), group.end(), 0);

      neighbor_query.add_group(group);

      RT_angles rt_angles(m_segments, angle_bound);
      rt_angles.add_group(group);
      rt_angles.make_bounds();

      SR_angles sr_angles(m_segments, neighbor_query, rt_angles);
      sr_angles.regularize();

      m_parallel_groups.clear();
      rt_angles.parallel_groups(
        std::back_inserter(m_parallel_groups));
    }

    void regularize_ordinates(const FT ordinate_bound) {
      std::cout << "Ordinate bound: " << ordinate_bound << std::endl;

      CGAL_assertion(m_parallel_groups.size() > 0);
      RT_ordinates rt_ordinates(m_segments, ordinate_bound);

      Neighbor_query neighbor_query(m_segments);
      for(const auto& parallel_group : m_parallel_groups) {
        neighbor_query.add_group(parallel_group);
        rt_ordinates.add_group(parallel_group);
      }

      SR_ordinates sr_ordinates(m_segments, neighbor_query, rt_ordinates);
      sr_ordinates.regularize();
    }

  private:
    Segment_range& m_segments;
    std::vector<Indices> m_parallel_groups;
  };

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_INTERNAL_REGULARIZATION_H