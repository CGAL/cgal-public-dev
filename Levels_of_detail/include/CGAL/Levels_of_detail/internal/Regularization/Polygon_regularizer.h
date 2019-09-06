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

// CGAL includes.
#include <CGAL/assertions.h>

// Internal includes.
#include <CGAL/Levels_of_detail/internal/utils.h>
#include <CGAL/Levels_of_detail/internal/struct.h>

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
    using Line_2 = typename Traits::Line_2;
    using Intersect_2 = typename Traits::Intersect_2;

    using FT_pair = std::pair<FT, FT>;
    using Size_pair = std::pair<std::size_t, std::size_t>;
    using Indices = std::vector<std::size_t>;

    Polygon_regularizer(
      const FT min_length,
      const FT angle_bound) :
    m_min_length(min_length),
    m_angle_bound(angle_bound),
    m_pi(static_cast<FT>(CGAL_PI)),
    m_angle_threshold(FT(5)) 
    { }

    void compute_longest_direction(
      const std::vector< std::vector<Segment_2> >& contours) {

      m_bounds.clear();
      m_bounds.resize(1);
      m_bounds[0] = std::make_pair(m_angle_bound, FT(90) - m_angle_bound);

      m_skip.clear();
      m_skip.resize(1);
      m_skip[0] = find_longest_segment(contours);

      m_longest.clear();
      m_longest.resize(1);
      m_longest[0] = contours[m_skip[0].first][m_skip[0].second];

      m_groups.clear();
      m_groups.resize(contours.size());
      for (std::size_t k = 0; k < contours.size(); ++k) {
        m_groups[k].resize(contours[k].size());
        for (std::size_t i = 0; i < contours[k].size(); ++i) {
          m_groups[k][i] = 0;
        }
      }
    }

    void regularize_contours(
      std::vector< std::vector<Segment_2> >& contours) {
      
      std::vector< std::vector<Segment_2> > tmp;
      for (std::size_t k = 0; k < contours.size(); ++k) {
        auto& contour = contours[k];

        rotate_contour(k, contour);
        correct_contour(k, contour);
        
        const bool success = connect_contour(contour);
        if (success)
          tmp.push_back(contour);
      }
      contours = tmp;
    }

    void rotate_contour(
      const std::size_t k,
      std::vector<Segment_2>& contour) {

      for (std::size_t i = 0; i < contour.size(); ++i) {
        const std::size_t gr_idx = m_groups[k][i];
        if (k == m_skip[gr_idx].first && i == m_skip[gr_idx].second) 
          continue;

        auto& si = contour[i];
        const auto& longest = m_longest[gr_idx];
        const auto& bounds = m_bounds[gr_idx];
        rotate_segment(longest, bounds, si);
      }
    }

    void correct_contour(
      const std::size_t k,
      std::vector<Segment_2>& contour) {

      const std::size_t n = contour.size();
      for (std::size_t i = 0; i < n; ++i) {
        const std::size_t gr_idx = m_groups[k][i];
        if (k == m_skip[gr_idx].first && i == m_skip[gr_idx].second) 
          continue;

        const std::size_t im = (i + n - 1) % n;
        const std::size_t ip = (i + 1) % n;

        auto& si = contour[i];
        const auto& sm = contour[im];
        const auto& sp = contour[ip];

        const FT length = internal::distance(si.source(), si.target());
        if (length < m_min_length)
          correct_segment(sm, si, sp);
      }
    }

    bool connect_contour(
      std::vector<Segment_2>& contour) {
      
      std::vector<Segment_2> initials, finals;
      remove_zero_length_segments(contour, initials);
      filter_out_wrong_segments(initials, finals);

      if (finals.size() < 4)
        return false;

      intersect_segments(finals);
      contour = finals;
      return true;
    }

    void remove_zero_length_segments(
      const std::vector<Segment_2>& contour,
      std::vector<Segment_2>& segments) {

      segments.clear();
      for (std::size_t i = 0; i < contour.size(); ++i) {
        const auto& si = contour[i];
        if (si.squared_length() > internal::tolerance<FT>())
          segments.push_back(si);
      }
    }

    void filter_out_wrong_segments(
      const std::vector<Segment_2>& contour,
      std::vector<Segment_2>& segments) {
      
      segments.clear();
      const std::size_t n = contour.size();
      const std::size_t start = find_initial_index(contour);
      std::size_t i = start;
      std::vector< std::vector<Segment_2> > tmp(1);
      do {

        tmp[0].clear(); std::size_t j = i;
        bool next_is_parallel = false;
        do {

          const std::size_t jm = (j + n - 1) % n;
          const std::size_t jp = (j + 1) % n;

          const auto& sj = contour[j];
          const auto& sm = contour[jm];
          const auto& sp = contour[jp];

          tmp[0].push_back(sj);
          const auto pair = check_segment(sm, sj, sp);
          next_is_parallel = pair.second;
          j = jp;

        } while (next_is_parallel);
        
        const auto data = find_longest_segment(tmp);
        segments.push_back(tmp[data.first][data.second]);
        i = j;

      } while (i != start);
    }

    std::size_t find_initial_index(
      const std::vector<Segment_2>& contour) {

      const std::size_t n = contour.size();
      for (std::size_t i = 0; i < n; ++i) {
        
        const std::size_t im = (i + n - 1) % n;
        const std::size_t ip = (i + 1) % n;
        
        const auto& si = contour[i];
        const auto& sm = contour[im];
        const auto& sp = contour[ip];

        const auto pair = check_segment(sm, si, sp);
        const bool previous_is_orthogonal = !(pair.first);
        if (previous_is_orthogonal) return i;
      }
    }

    void intersect_segments(
      std::vector<Segment_2>& segments) {

      const std::size_t n = segments.size();
      for (std::size_t i = 0; i < n; ++i) {
        
        const std::size_t im = (i + n - 1) % n;
        const std::size_t ip = (i + 1) % n;
        
        auto& si = segments[i];
        const auto& sm = segments[im];
        const auto& sp = segments[ip];
        
        intersect_segment(sm, si, sp);
      }
    }

  private:
    const FT m_min_length;
    const FT m_angle_bound;

    const FT m_pi;
    const FT m_angle_threshold;

    std::vector<FT_pair> m_bounds;
    std::vector<Size_pair> m_skip;
    std::vector<Segment_2> m_longest;
    std::vector<Indices> m_groups;

    Size_pair find_longest_segment(
      const std::vector< std::vector<Segment_2> >& contours) {

      std::size_t con_idx = std::size_t(-1);
      std::size_t seg_idx = std::size_t(-1);

      FT max_length = -FT(1);
      for (std::size_t k = 0; k < contours.size(); ++k) {
        for (std::size_t i = 0; i < contours[k].size(); ++i) {
          
          const FT length = contours[k][i].squared_length();
          if (length > max_length) {

            max_length = length;
            con_idx = k;
            seg_idx = i;
          }
        }
      }
      return std::make_pair(con_idx, seg_idx);
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

    std::pair<bool, bool> check_segment(
      const Segment_2& sm, const Segment_2& si, const Segment_2& sp) {

      // Source.
      const FT angle_mi = angle_degree_2(sm, si);
      const FT angle_mi_2 = get_angle_2(angle_mi);

      // Target.
      const FT angle_pi = angle_degree_2(si, sp);
      const FT angle_pi_2 = get_angle_2(angle_pi);

      const bool source_cond = ( CGAL::abs(angle_mi_2) <= m_angle_threshold );
      const bool target_cond = ( CGAL::abs(angle_pi_2) <= m_angle_threshold );

      return std::make_pair(source_cond, target_cond);
    }

    void intersect_segment(
      const Segment_2& sm, Segment_2& si, const Segment_2& sp) {

      Point_2 source = si.source();
      Point_2 target = si.target();

      const Line_2 line_1 = Line_2(sm.source(), sm.target());
      const Line_2 line_2 = Line_2(si.source(), si.target());
      const Line_2 line_3 = Line_2(sp.source(), sp.target());

      /* // use in case projection does not work!
      const bool success1 = intersect_2(line_1, line_2, source);
      const bool success2 = intersect_2(line_2, line_3, target);

      if (!success1) source = si.source();
      if (!success2) target = si.target(); */

      source = line_1.projection(si.source());
      target = line_3.projection(si.target());

      si = Segment_2(source, target);
    }

    bool intersect_2(
      const Line_2& line_1, const Line_2& line_2,
      Point_2& in_point) {
      
      typename std::result_of<Intersect_2(Line_2, Line_2)>::type result 
      = CGAL::intersection(line_1, line_2);
      if (result) {
        if (const Line_2* line = boost::get<Line_2>(&*result)) 
          return false;
        else {
          const Point_2* point = boost::get<Point_2>(&*result);
          in_point = *point; return true;
        }
      }
      return false;
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
  };

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_INTERNAL_POLYGON_REGULARIZER_H
