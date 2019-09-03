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

#ifndef CGAL_LEVELS_OF_DETAIL_INTERNAL_POLYGON_REGULARIZER_H
#define CGAL_LEVELS_OF_DETAIL_INTERNAL_POLYGON_REGULARIZER_H

#include <CGAL/license/Levels_of_detail.h>

// STL includes.
#include <vector>
#include <utility>

// Internal includes.
#include <CGAL/Levels_of_detail/internal/utils.h>
#include <CGAL/Levels_of_detail/internal/struct.h>

// Testing.
#include "../../../../../test/Levels_of_detail/include/Saver.h"

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

  template<typename GeomTraits>
  class Polygon_regularizer {

  public:
    using Traits = GeomTraits;

    using FT = typename Traits::FT;
    using Point_2 = typename Traits::Point_2;
    using Segment_2 = typename Traits::Segment_2;
    using Vector_2 = typename Traits::Vector_2;

    using Saver = Saver<Traits>;
    using Color = CGAL::Color;

    Polygon_regularizer(
      const FT min_length,
      const FT angle_bound) :
    m_min_length(min_length),
    m_angle_bound(angle_bound),
    m_pi(static_cast<FT>(CGAL_PI)),
    m_bound_1(m_pi / FT(16)),
    m_bound_2(m_pi / FT(2) - m_bound_1)
    { }

    void regularize(std::vector<Segment_2>& contour) {

      const std::size_t n = contour.size();
      const std::size_t ref_idx = find_longest_segment(contour);
      const auto& rs = contour[ref_idx];
      
      for (std::size_t i = 0; i < n; ++i) {
        
        const std::size_t im = (i + n - 1) % n;
        const std::size_t im = (i + 1) % n;
        
        const auto& si = contour[i];
        const auto& sm = contour[im];
        const auto& sp = contour[ip];

        const FT length = internal::distance(si.source(), si.target());
        if (length > m_min_length)
          handle_long_segment(rs, sm, si, sp);
        else
          handle_short_segment(rs, sm, si, sp);
      }
    }

  private:
    const FT m_min_length;
    const FT m_angle_bound;

    const FT m_pi;
    const FT m_bound_1;
    const FT m_bound_2;

    std::size_t find_longest_segment(
      const std::vector<Segment_2>& contour) {

      std::size_t idx = std::size_t(-1);
      FT max_length = -FT(1);
      for (std::size_t i = 0; i < contour.size(); ++i) {
        const FT length = contour[i].squared_length();
        if (length > max_length) {
          max_length = length;
          idx = i;
        }
      }
      return idx;
    }

    void handle_long_segment(
      const Segment_2& longest,
      Segment_2& sm, Segment_2& si, Segment_2& sp) {

      const Vector_2 v1 = si.to_vector();
      const Vector_2 v2 = longest.to_vector();

      FT angle_rad;
      internal::compute_angle_2(v1, v2, angle_rad);
      // const FT angle_deg = angle_rad * FT(180) / m_pi;

      if (angle_rad <= m_bound_1)
        rotate(angle_rad, sm, si, sp);
      if (angle_rad >= m_bound_2)
        rotate(m_pi - angle_rad, sm, si, sp);
    }

    void handle_short_segment(
      const Segment_2& longest,
      Segment_2& sm, Segment_2& si, Segment_2& sp) {

      
    }

    void rotate(
      const FT angle_2,
      Segment_2& sm, Segment_2& si, Segment_2& sp) {

      Point_2 source_i = si.source;
      Point_2 target_i = si.target;
      const Point_2 b = internal::middle_point(source_i, target_i);
      internal::rotate_point_2(angle_2, b, source_i);
      internal::rotate_point_2(angle_2, b, target_i);

      sm = Segment_2(sm.source(), source_i);
      si = Segment_2(source_i, target_i);
      sp = Segment_2(target_i, sp.target());
    }
  };

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_INTERNAL_POLYGON_REGULARIZER_H
