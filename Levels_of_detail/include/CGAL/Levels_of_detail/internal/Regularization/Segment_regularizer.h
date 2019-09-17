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

#ifndef CGAL_LEVELS_OF_DETAIL_INTERNAL_SEGMENT_REGULARIZER_H
#define CGAL_LEVELS_OF_DETAIL_INTERNAL_SEGMENT_REGULARIZER_H

#include <CGAL/license/Levels_of_detail.h>

// STL includes.
#include <vector>
#include <utility>
#include <algorithm>

// CGAL includes.
#include <CGAL/assertions.h>

// Internal includes.
#include <CGAL/Levels_of_detail/internal/utils.h>
#include <CGAL/Levels_of_detail/internal/struct.h>

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

  template<typename GeomTraits>
  class Segment_regularizer {

  public:
    using Traits = GeomTraits;

    using FT = typename Traits::FT;
    using Point_2 = typename Traits::Point_2;
    using Segment_2 = typename Traits::Segment_2;
    using Vector_2 = typename Traits::Vector_2;
    using Line_2 = typename Traits::Line_2;

    using FT_pair = std::pair<FT, FT>;
    using Indices = std::vector<std::size_t>;

    Segment_regularizer(
      const FT min_length,
      const FT angle_bound,
      const FT ordinate_bound) :
    m_min_length(min_length),
    m_angle_bound(angle_bound),
    m_ordinate_bound(ordinate_bound),
    m_pi(static_cast<FT>(CGAL_PI)),
    m_angle_threshold(FT(5))
    { }

    void compute_longest_direction(
      const std::vector<Segment_2>& outer_segments,
      const std::vector< std::vector<Segment_2> >& contours) {

      m_bounds.clear();
      m_bounds.resize(1);
      m_bounds[0] = std::make_pair(m_angle_bound, FT(90) - m_angle_bound);

      const std::size_t seg_idx = find_longest_segment(outer_segments);

      m_longest.clear();
      m_longest.resize(1);
      m_longest[0] = outer_segments[seg_idx];

      m_groups.clear();
      m_groups.resize(contours.size());
      for (std::size_t k = 0; k < contours.size(); ++k) {
        m_groups[k].resize(contours[k].size(), std::size_t(-1));
        for (std::size_t i = 0; i < contours[k].size(); ++i) {
          m_groups[k][i] = 0;
        }
      }
    }

    void regularize_contours(
      std::vector< std::vector<Segment_2> >& contours) {
      
      if (m_min_length == FT(0))
        return;

      for (std::size_t k = 0; k < contours.size(); ++k) {
        auto& contour = contours[k];
        rotate_contour(k, contour);

        /* if (contour.size() == 2)
          correct_contour_2(k, contour); */

        if (contour.size() >= 3)
          correct_contour_n(k, contour);
      }
    }

    void merge_closest(std::vector<Segment_2>& segments) {
      
      std::vector<Segment_2> merged;

      std::vector<bool> states(segments.size(), false);
      std::vector<Segment_2> group; std::size_t num_groups = 0;
      for (std::size_t i = 0; i < segments.size(); ++i) {
        if (states[i]) continue;
        
        group.clear();
        const auto& segment = segments[i];
        states[i] = true;
        group.push_back(segment);
        const auto& p = segment.source();

        for (std::size_t j = 0; j < segments.size(); ++j) {
          if (states[j]) continue;
          
          const FT angle = angle_degree_2(segment, segments[j]);
          const FT angle_2 = get_angle_2(angle);

          if (CGAL::abs(angle_2) <= m_angle_threshold) {
            const Line_2 line = Line_2(
              segments[j].source(), segments[j].target());
            const auto q = line.projection(p);
            const FT dist = internal::distance(p, q);
            if (dist <= m_ordinate_bound) {
              group.push_back(segments[j]);
              states[j] = true;
            }
          }
        }

        auto central = find_central_segment(group);
        const Line_2 line = Line_2(central.source(), central.target());
        
        FT sum_length = FT(0);
        std::vector<Point_2> points;
        for (const auto& seg : group) {
          const Point_2 p = line.projection(seg.source());
          const Point_2 q = line.projection(seg.target());
          points.push_back(p);
          points.push_back(q);
          sum_length += internal::distance(p, q);
        }
        update_segment(points, central);
        merged.push_back(central);

        ++num_groups;
      }
      segments = merged;
      std::cout << "Num collinear groups: " << num_groups << std::endl;
    }

  private:
    const FT m_min_length;
    const FT m_angle_bound;
    const FT m_ordinate_bound;

    const FT m_pi;
    const FT m_angle_threshold;

    std::vector<FT_pair> m_bounds;
    std::vector<Segment_2> m_longest;
    std::vector<Indices> m_groups;

    std::size_t find_longest_segment(
      const std::vector<Segment_2>& segments) {

      std::size_t seg_idx = std::size_t(-1);
      FT max_length = -FT(1);
      for (std::size_t i = 0; i < segments.size(); ++i) {
          
        const FT length = segments[i].squared_length();
        if (length > max_length) {

          max_length = length;
          seg_idx = i;
        }
      }
      return seg_idx;
    }

    void rotate_contour(
      const std::size_t k,
      std::vector<Segment_2>& contour) {

      for (std::size_t i = 0; i < contour.size(); ++i) {
        const std::size_t gr_idx = m_groups[k][i];
        if (gr_idx == std::size_t(-1))
          continue;

        auto& si = contour[i];
        const auto& longest = m_longest[gr_idx];
        const auto& bounds = m_bounds[gr_idx];
        rotate_segment(longest, bounds, si);
      }
    }

    void rotate_segment(
      const Segment_2& longest,
      const FT_pair& bounds,
      Segment_2& si) {

      const FT angle = angle_degree_2(longest, si);
      const FT angle_2 = get_angle_2(angle);

      if (CGAL::abs(angle_2) <= bounds.first) {
        rotate(angle, FT(180), longest, si); // parallel case
        return;
      }
      if (CGAL::abs(angle_2) >= bounds.second) {
        rotate(angle, FT(90), longest, si); // orthogonal case
        return;
      }
    }

    FT angle_degree_2(
      const Segment_2& longest, const Segment_2& si) {

      const Vector_2 v1 = si.to_vector();
      const Vector_2 v2 = -longest.to_vector();

		  const FT det = CGAL::determinant(v1, v2);
		  const FT dot = CGAL::scalar_product(v1, v2);
      const FT angle_rad = static_cast<FT>(
        std::atan2(CGAL::to_double(det), CGAL::to_double(dot)));
      const FT angle_deg = angle_rad * FT(180) / m_pi; 
      return angle_deg;
    }

    FT get_angle_2(const FT angle) {
      
      FT angle_2 = angle;
      if (angle_2 > FT(90)) angle_2 = FT(180) - angle_2;
      else if (angle_2 < -FT(90)) angle_2 = FT(180) + angle_2;
      return angle_2;
    }

    void rotate(
      const FT angle_2, 
      const FT ref_angle_2,
      const Segment_2& longest, 
      Segment_2& si) {

      FT angle = angle_2;
      if (angle < FT(0)) angle = angle + ref_angle_2;
      else if (angle > FT(0)) angle = angle - ref_angle_2;

      Point_2 source_i = si.source();
      Point_2 target_i = si.target();
      const Point_2 b = internal::middle_point_2(source_i, target_i);
      const FT angle_rad = angle * m_pi / FT(180); 
      internal::rotate_point_2(angle_rad, b, source_i);
      internal::rotate_point_2(angle_rad, b, target_i);
      si = Segment_2(source_i, target_i);
    }

    void correct_contour_2(
      const std::size_t k,
      std::vector<Segment_2>& contour) {

      const std::size_t seg_idx = find_longest_segment(contour);
      
      const std::size_t i = (seg_idx == 0 ? 1 : 0);
      const std::size_t j =  seg_idx;
      
      const std::size_t gr_idx = m_groups[k][i];
      if (gr_idx == std::size_t(-1))
        return;

      auto& si = contour[i];
      const auto& sj = contour[j];

      FT angle_ij = FT(0);
      if (i > j) angle_ij = angle_degree_2(sj, si);
      else angle_ij = angle_degree_2(si, sj);

      const FT angle_ij_2 = get_angle_2(angle_ij);
      if (CGAL::abs(angle_ij_2) <= m_angle_threshold) {
        
        FT angle = FT(0);
        if (i > j) angle = angle_degree_2(sj, si);
        else angle = angle_degree_2(si, sj);
        
        rotate(angle, FT(90), sj, si); // orthogonal case
        return;
      }
    }

    void correct_contour_n(
      const std::size_t k,
      std::vector<Segment_2>& contour) {

      const std::size_t n = contour.size();
      for (std::size_t i = 1; i < n - 1; ++i) {
        const std::size_t gr_idx = m_groups[k][i];
        if (gr_idx == std::size_t(-1))
          continue;

        const std::size_t im = i - 1;
        const std::size_t ip = i + 1;

        auto& si = contour[i];
        const auto& sm = contour[im];
        const auto& sp = contour[ip];

        const FT length = internal::distance(si.source(), si.target());
        if (length <= m_min_length)
          correct_segment(sm, si, sp);
      }
    }

    void correct_segment(
      const Segment_2& sm, Segment_2& si, const Segment_2& sp) {

      const FT angle_mp = angle_degree_2(sm, sp);
      const FT angle_mp_2 = get_angle_2(angle_mp);

      if (CGAL::abs(angle_mp_2) <= m_angle_threshold) {
        const FT angle = angle_degree_2(sm, si);
        rotate(angle, FT(90), sm, si); // orthogonal case
        return;
      }
    }

    Segment_2 find_central_segment(
      const std::vector<Segment_2>& segments) {

      Point_2 source, target;
      FT x1 = FT(0), y1 = FT(0);
      FT x2 = FT(0), y2 = FT(0);
      for (const auto& segment : segments) {
        x1 += segment.source().x();
        x2 += segment.target().x();

        y1 += segment.source().y();
        y2 += segment.target().y();
      }

      const FT size = static_cast<FT>(segments.size());
      x1 /= size; y1 /= size;
      x2 /= size; y2 /= size;

      return Segment_2(Point_2(x1, y1), Point_2(x2, y2));
    }

    void update_segment(
      const std::vector<Point_2>& points,
      Segment_2& segment) {

      FT min_proj_value =  internal::max_value<FT>();
      FT max_proj_value = -internal::max_value<FT>();

      const Vector_2 ref_vector = segment.to_vector();
      Point_2 ref_point;
      internal::compute_barycenter_2(points, ref_point);
      
      Point_2 p, q;
      for (const auto& point : points) {
        const Vector_2 curr_vector(ref_point, point);
        const FT value = CGAL::scalar_product(curr_vector, ref_vector);
        
        if (value < min_proj_value) {
          min_proj_value = value;
          p = point; }
        if (value > max_proj_value) {
          max_proj_value = value;
          q = point; }
      }
      segment = Segment_2(p, q);
    }
  };

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_INTERNAL_SEGMENT_REGULARIZER_H
