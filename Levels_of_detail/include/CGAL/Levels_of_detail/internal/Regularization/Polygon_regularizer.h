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
#include <algorithm>

// CGAL includes.
#include <CGAL/assertions.h>

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
    using Point_3 = typename Traits::Point_3;
    using Segment_2 = typename Traits::Segment_2;
    using Vector_2 = typename Traits::Vector_2;
    using Line_2 = typename Traits::Line_2;
    using Intersect_2 = typename Traits::Intersect_2;

    using FT_pair = std::pair<FT, FT>;
    using Size_pair = std::pair<std::size_t, std::size_t>;
    using Indices = std::vector<std::size_t>;
    using Seg_pair = std::pair<Segment_2, bool>;

    Polygon_regularizer(
      const FT min_length,
      const FT angle_bound,
      const FT ordinate_bound) :
    m_min_length(min_length),
    m_angle_bound(angle_bound),
    m_ordinate_bound(ordinate_bound * FT(2)),
    m_pi(static_cast<FT>(CGAL_PI)),
    m_angle_threshold(FT(5)),
    m_bound_min(m_angle_bound),
    m_bound_max(FT(90) - m_bound_min)
    { }

    void compute_multiple_directions(
      const std::vector< std::vector<Segment_2> >& input_contours) {

      std::vector< std::vector<Seg_pair> > contours;
      create_internal_contours(input_contours, contours);
      get_multiple_directions(contours);

      if (m_longest.size() != 0) {

        unify_along_contours(contours);
        correct_directions(contours);
        readjust_directions(contours);
      }

      if (m_longest.size() == 0)
        compute_longest_direction(contours);

      std::cout << "Num outer directions: " << m_longest.size() << std::endl;
    }

    void create_internal_contours(
      const std::vector< std::vector<Segment_2> >& input,
      std::vector< std::vector<Seg_pair> >& output) {

      output.clear();
      std::vector<Seg_pair> segments;
      for (const auto& contour : input) {

        segments.clear();
        for (const auto& segment : contour) {
          const auto& s = segment.source();
          const auto& t = segment.target();

          if (internal::distance(s, t) >= m_min_length * FT(2))
            segments.push_back(std::make_pair(segment, true));
          else
            segments.push_back(std::make_pair(segment, false));
        }
        output.push_back(segments);
      }
    }

    void unify_along_contours(
      const std::vector< std::vector<Seg_pair> >& contours) {

      for (std::size_t k = 0; k < contours.size(); ++k) {
        for (std::size_t i = 0; i < contours[k].size(); ++i) {
          if (contours[k][i].second)
            continue;

          if (m_groups[k][i] == std::size_t(-1)) {

            const std::size_t m = contours[k].size();
            std::size_t im = (i + m - 1) % m;
            std::size_t ip = (i + 1) % m;

            bool stop = false;
            std::size_t max_count = 0;
            do {

              if (contours[k][im].second) {
                m_groups[k][i] = m_groups[k][im]; break;
              }

              if (contours[k][ip].second) {
                m_groups[k][i] = m_groups[k][ip]; break;
              }

              im = (im + m - 1) % m;
              ip = (ip + 1) % m;

              if (im == i || ip == i) stop = true;
              ++max_count;

            } while (!stop && max_count < m * 2);
            if (stop || max_count > m * 2)
              m_groups[k][i] = 0;
          }
        }
      }
    }

    void unify_along_contours(
      const std::vector<Seg_pair>& contour,
      Indices& group) {

      const std::size_t m = contour.size();
      for (std::size_t i = 0; i < m; ++i) {

        if (contour[i].second) continue;
        if (group[i] == std::size_t(-1)) {

          std::size_t im = std::size_t(-1);
          if (i > 0 && i <= m) im = i - 1;
          else im = std::size_t(-1);

          std::size_t ip = std::size_t(-1);
          if (i < m - 1 && i >= 0) ip = i + 1;
          else ip = std::size_t(-1);

          bool stop = false;
          std::size_t max_count = 0;
          do {

            if (im != std::size_t(-1)) {
              if (contour[im].second) {
                if (group[im] == std::size_t(-1)) {
                  group[i] = 0;
                } else {
                  group[i] = group[im];
                }
                break;
              }
            }

            if (ip != std::size_t(-1)) {
              if (contour[ip].second) {
                if (group[ip] == std::size_t(-1)) {
                  group[i] = 0;
                } else {
                  group[i] = group[ip];
                }
                break;
              }
            }

            if (im != std::size_t(-1)) {
              std::size_t tmp = im;
              if (tmp > 0 && tmp < m) im = tmp - 1;
              else im = std::size_t(-1);
            }

            if (ip != std::size_t(-1)) {
              std::size_t tmp = ip;
              if (tmp < m - 1 && tmp >= 0) ip = tmp + 1;
              else ip = std::size_t(-1);
            }

            if (im == std::size_t(-1) && ip == std::size_t(-1)) stop = true;
            ++max_count;

          } while (!stop && max_count < m);
          if (stop || max_count > m)
            group[i] = 0;
        }
      }
    }

    void correct_directions(
      const std::vector< std::vector<Seg_pair> >& contours) {

      for (std::size_t k = 0; k < contours.size(); ++k) {
        const std::size_t n = contours[k].size();

        Indices group; group.reserve(n);
        for (std::size_t i = 0; i < n; ++i) {

          const std::size_t im = (i + n - 1) % n;
          const std::size_t ip = (i + 1) % n;

          const std::size_t gm = m_groups[k][im];
          const std::size_t gi = m_groups[k][i];
          const std::size_t gp = m_groups[k][ip];

          const bool skipped = ( k == m_skip[gi].first && i == m_skip[gi].second );
          if (gm == gp && gi != gm && !skipped)
            group.push_back(gm);
          else
            group.push_back(gi);
        }
        m_groups[k] = group;
      }
    }

    void correct_directions(
      const std::vector<Seg_pair>& contour,
      Indices& group) {

      const std::size_t n = contour.size();
      if (n <= 2) return;

      Indices clean; clean.reserve(n);
      clean.push_back(group[0]);
      for (std::size_t i = 1; i < n - 1; ++i) {

        const std::size_t im = i - 1;
        const std::size_t ip = i + 1;

        const std::size_t gm = group[im];
        const std::size_t gi = group[i];
        const std::size_t gp = group[ip];

        if (gm == gp && gi != gm)
          clean.push_back(gm);
        else
          clean.push_back(gi);
      }
      clean.push_back(group[n - 1]);
      group = clean;
    }

    void readjust_directions(
      const std::vector< std::vector<Seg_pair> >& contours) {

      std::vector<FT> angles, counts;
      create_average_angles(contours, angles, counts);

      Segment_2 stub;
      for (std::size_t k = 0; k < angles.size(); ++k) {
        angles[k] /= counts[k];
        const FT angle = angles[k];
        rotate(angle, FT(0), stub, m_longest[k]);
      }
    }

    template<typename T>
    void make_default_groups(
      const std::vector< std::vector<T> >& contours,
      const std::size_t value,
      std::vector<Indices>& groups) {

      groups.clear(); groups.resize(contours.size());
      for (std::size_t k = 0; k < contours.size(); ++k)
        groups[k].resize(contours[k].size(), value);
    }

    FT angle_degree_2(
      const Segment_2& longest, const Segment_2& segment) {

      const Vector_2 v1 =  segment.to_vector();
      const Vector_2 v2 = -longest.to_vector();

		  const FT det = CGAL::determinant(v1, v2);
		  const FT dot = CGAL::scalar_product(v1, v2);
      const FT angle_rad = static_cast<FT>(
        std::atan2(CGAL::to_double(det), CGAL::to_double(dot)));
      const FT angle_deg = angle_rad * FT(180) / m_pi;
      return angle_deg;
    }

    const FT get_bound_min() const {
      return m_bound_min;
    }

    const FT get_bound_max() const {
      return m_bound_max;
    }

    FT get_angle_2(const FT angle) {

      FT angle_2 = angle;
      if (angle_2 > FT(90)) angle_2 = FT(180) - angle_2;
      else if (angle_2 < -FT(90)) angle_2 = FT(180) + angle_2;
      return angle_2;
    }

    void set_data(
      const std::vector<FT_pair>& bounds,
      const std::vector<Size_pair>& skip,
      const std::vector<Segment_2>& longest,
      const std::vector<Indices>& groups) {

      m_bounds = bounds;
      m_skip = skip;
      m_longest = longest;
      m_groups = groups;
    }

    const std::vector<Segment_2>& get_longest() const {
      return m_longest;
    }

    const std::vector<Indices>& get_groups() const {
      return m_groups;
    }

    void regularize_contours(
      std::vector< std::vector<Segment_2> >& contours) {

      if (m_angle_bound == FT(0))
        return;

      std::vector< std::vector<Segment_2> > initials, finals;
      initials.reserve(contours.size());
      for (const auto& contour : contours)
        initials.push_back(contour);
      m_saved = initials;

      for (std::size_t k = 0; k < initials.size(); ++k) {
        auto contour = initials[k];

        rotate_contour(k, contour);

        /* correct_contour(k, contour); */

        optimize_contour(contour);

        const bool success = connect_contour(contour);
        if (success)
          finals.push_back(contour);
      }
      contours = finals;
    }

    void regularize_polyline(
      std::vector<Segment_2>& contour) {

      if (m_angle_bound == FT(0))
        return;

      auto init = contour;
      rotate_polyline_contour(init);
      optimize_contour(init);
      const bool success = connect_polyline(init);
      if (success)
        contour = init;
    }

    void get_directions(
      std::vector<Segment_2>& directions) {

      directions.clear();
      directions = m_longest;
    }

  private:
    const FT m_min_length;
    const FT m_angle_bound;
    const FT m_ordinate_bound;

    const FT m_pi;
    const FT m_angle_threshold;
    const FT m_bound_min, m_bound_max;

    std::vector<FT_pair> m_bounds;
    std::vector<Size_pair> m_skip;
    std::vector<Segment_2> m_longest;
    std::vector<Indices> m_groups;

    std::vector< std::vector<Segment_2> > m_saved;

    void compute_longest_direction(
      const std::vector< std::vector<Seg_pair> >& contours) {

      m_bounds.clear(); m_bounds.resize(1);
      m_bounds[0] = std::make_pair(FT(45), FT(45));

      const auto longest_pair = find_longest_segment(contours);

      m_skip.clear(); m_skip.resize(1);
      m_skip[0] = longest_pair;

      m_longest.clear(); m_longest.resize(1);
      m_longest[0] = contours[longest_pair.first][longest_pair.second].first;

      make_default_groups(contours, 0, m_groups);
    }

    void get_multiple_directions(
      const std::vector< std::vector<Seg_pair> >& contours) {

      std::vector<Size_pair> input;
      for (std::size_t k = 0; k < contours.size(); ++k)
        for (std::size_t i = 0; i < contours[k].size(); ++i)
          input.push_back(std::make_pair(k, i));

      sort_input(contours, input);
      std::vector<bool> states(input.size(), false);

      m_bounds.clear(); m_skip.clear();
      m_longest.clear(); m_groups.clear();
      make_default_groups(contours, std::size_t(-1), m_groups);

      bool apply = true; std::size_t gr_idx = 0;
      do {
        apply = get_next_direction(
          contours, input, gr_idx, states);
        ++gr_idx;
      } while (apply);
    }

    void sort_input(
      const std::vector< std::vector<Seg_pair> >& contours,
      std::vector<Size_pair>& input) {

      std::sort(input.begin(), input.end(),
      [&contours](const Size_pair& a, const Size_pair& b) -> bool {
        const FT length_1 = (contours[a.first][a.second]).first.squared_length();
        const FT length_2 = (contours[b.first][b.second]).first.squared_length();
        return length_1 > length_2;
      });
    }

    bool get_next_direction(
      const std::vector< std::vector<Seg_pair> >& contours,
      const std::vector<Size_pair>& input,
      const std::size_t gr_idx,
      std::vector<bool>& states) {

      std::size_t longest_idx = std::size_t(-1);
      for (std::size_t i = 0; i < states.size(); ++i) {
        if (!states[i] && contours[input[i].first][input[i].second].second) {
          longest_idx = i; break;
        }
      }
      if (longest_idx == std::size_t(-1))
        return false;

      const auto& longest_pair = input[longest_idx];
      const Segment_2& longest_segment =
        contours[longest_pair.first][longest_pair.second].first;

      for (std::size_t i = 0; i < states.size(); ++i) {
        if (i == longest_idx) {
          m_groups[longest_pair.first][longest_pair.second] = gr_idx;
          states[i] = true; continue;
        }

        const auto& pair = input[i];
        if (!states[i] && contours[pair.first][pair.second].second) {
          const auto& segment = contours[pair.first][pair.second].first;

          const FT angle = angle_degree_2(longest_segment, segment);
          const FT angle_2 = get_angle_2(angle);

          if (
            (CGAL::abs(angle_2) <= m_bound_min) ||
            (CGAL::abs(angle_2) >= m_bound_max) )  {

            m_groups[pair.first][pair.second] = gr_idx;
            states[i] = true; continue;
          }
        }
      }

      m_longest.push_back(longest_segment);
      m_bounds.push_back(std::make_pair(FT(45), FT(45)));
      m_skip.push_back(longest_pair);

      return true;
    }

    Size_pair find_longest_segment(
      const std::vector< std::vector<Seg_pair> >& contours) {

      std::size_t con_idx = std::size_t(-1);
      std::size_t seg_idx = std::size_t(-1);

      FT max_length = -FT(1);
      for (std::size_t k = 0; k < contours.size(); ++k) {
        for (std::size_t i = 0; i < contours[k].size(); ++i) {

          const auto& segment = contours[k][i].first;
          const FT length = segment.squared_length();
          if (length > max_length) {

            max_length = length;
            con_idx = k; seg_idx = i;
          }
        }
      }
      return std::make_pair(con_idx, seg_idx);
    }

    Segment_2 find_longest_segment(
      const std::vector<Segment_2>& segments) {

      FT max_length = -FT(1);
      std::size_t seg_idx = std::size_t(-1);

      for (std::size_t i = 0; i < segments.size(); ++i) {
        const auto& segment = segments[i];
        const FT length = segment.squared_length();
        if (length > max_length) {

          max_length = length;
          seg_idx = i;
        }
      }
      return segments[seg_idx];
    }

    void create_average_angles(
      const std::vector< std::vector<Seg_pair> >& contours,
      std::vector<FT>& angles,
      std::vector<FT>& counts) {

      angles.clear();
      angles.resize(m_longest.size(), FT(0));

      counts.clear();
      counts.resize(m_longest.size(), FT(0));

      for (std::size_t k = 0; k < m_groups.size(); ++k) {
        for (std::size_t i = 0; i < m_groups[k].size(); ++i) {

          if (!contours[k][i].second) continue;
          const std::size_t gr_idx = m_groups[k][i];

          const auto& si = m_longest[gr_idx];
          const auto& sj = contours[k][i].first;

          const auto di = internal::compute_direction(si);
          const auto dj = internal::compute_direction(sj);

          const FT oi = internal::compute_orientation(di);
          const FT oj = internal::compute_orientation(dj);

          const FT mes_ij = oi - oj;
          const double mes90 = std::floor(CGAL::to_double(mes_ij / FT(90)));

          const FT to_lower = FT(90) *  static_cast<FT>(mes90)          - mes_ij;
          const FT to_upper = FT(90) * (static_cast<FT>(mes90) + FT(1)) - mes_ij;

          const FT angle = CGAL::abs(to_lower) < CGAL::abs(to_upper) ? to_lower : to_upper;

          angles[gr_idx] += angle;
          counts[gr_idx] += FT(1);
        }
      }
    }

    void rotate_contour(
      const std::size_t k,
      std::vector<Segment_2>& contour) {

      for (std::size_t i = 0; i < contour.size(); ++i) {
        const std::size_t gr_idx = m_groups[k][i];

        if (gr_idx == std::size_t(-1))
          continue;
        if (k == m_skip[gr_idx].first && i == m_skip[gr_idx].second)
          continue;

        auto& segment = contour[i];
        const auto& longest_segment = m_longest[gr_idx];
        const auto& bounds = m_bounds[gr_idx];

        const bool success = rotate_segment(longest_segment, bounds, segment);
        if (!success)
          m_groups[k][i] = std::size_t(-1);
      }
    }

    void rotate_polyline_contour(
      std::vector<Segment_2>& contour) {

      for (std::size_t i = 0; i < contour.size(); ++i) {
        const std::size_t gr_idx = m_groups[0][i];

        auto& segment = contour[i];
        if (gr_idx == std::size_t(-1)) {
          const auto& longest_segment = m_longest[0];
          const auto& bounds = m_bounds[0];

          const FT angle = angle_degree_2(longest_segment, segment);
          const FT angle_2 = get_angle_2(angle);
          rotate(angle, FT(180), longest_segment, segment);

        } else {
          const auto& longest_segment = m_longest[gr_idx];
          const auto& bounds = m_bounds[gr_idx];
          rotate_segment(longest_segment, bounds, segment);
        }
      }
    }

    bool rotate_segment(
      const Segment_2& longest_segment,
      const FT_pair& bounds,
      Segment_2& segment) {

      const FT angle = angle_degree_2(longest_segment, segment);
      const FT angle_2 = get_angle_2(angle);

      if (CGAL::abs(angle_2) <= bounds.first)
        rotate(angle, FT(180), longest_segment, segment); // parallel case
      if (CGAL::abs(angle_2) >= bounds.second)
        rotate(angle, FT(90), longest_segment, segment); // orthogonal case
      return true;
    }

    void rotate(
      const FT angle_2,
      const FT ref_angle_2,
      const Segment_2& longest_segment,
      Segment_2& segment) {

      FT angle = angle_2;
      if (angle < FT(0)) angle = angle + ref_angle_2;
      else if (angle > FT(0)) angle = angle - ref_angle_2;

      Point_2 source = segment.source();
      Point_2 target = segment.target();

      const Point_2 b = internal::middle_point_2(source, target);
      const FT angle_rad = angle * m_pi / FT(180);

      internal::rotate_point_2(angle_rad, b, source);
      internal::rotate_point_2(angle_rad, b, target);

      segment = Segment_2(source, target);
    }

    void correct_contour(
      const std::size_t k,
      std::vector<Segment_2>& contour) {

      const std::size_t n = contour.size();
      std::vector<Segment_2> segments;

      segments.reserve(n);
      for (std::size_t i = 0; i < n; ++i) {

        const auto& si = contour[i];
        const std::size_t gr_idx = m_groups[k][i];

        if (gr_idx == std::size_t(-1)) {
          segments.push_back(si); continue;
        }

        if (k == m_skip[gr_idx].first && i == m_skip[gr_idx].second) {
          segments.push_back(si); continue;
        }

        const std::size_t im = (i + n - 1) % n;
        const std::size_t ip = (i + 1) % n;

        auto ss = si;
        const auto& sm = contour[im];
        const auto& sp = contour[ip];

        const auto& inm = m_saved[k][im];
        const auto& inp = m_saved[k][ip];

        correct_segment(sm, ss, sp, inm, inp);
        segments.push_back(ss);
      }
      contour = segments;
    }

    void correct_segment(
      const Segment_2& sm, Segment_2& si, const Segment_2& sp,
      const Segment_2& inm, const Segment_2& inp) {

      const FT angle_mp   = angle_degree_2(sm, sp);
      const FT angle_mp_2 = get_angle_2(angle_mp);

      const FT angle_mi   = angle_degree_2(sm, si);
      const FT angle_mi_2 = get_angle_2(angle_mi);

      const FT angle_in   = angle_degree_2(inm, inp);
      const FT angle_in_2 = get_angle_2(angle_in);

      const FT length = internal::distance(si.source(), si.target());
      if (
        CGAL::abs(angle_mp_2) <= m_angle_threshold &&
        CGAL::abs(angle_mi_2) <= m_angle_threshold &&
        length <= m_min_length) {

        rotate(angle_mi, FT(90), sm, si); // orthogonal case
        return;
      }

      if (
        CGAL::abs(angle_mp_2) <= m_angle_threshold &&
        CGAL::abs(angle_mi_2) <= m_angle_threshold &&
        CGAL::abs(angle_in_2) <= m_angle_threshold &&
        length > m_min_length) {

        rotate(angle_mi, FT(90), sm, si); // orthogonal case
        return;
      }
    }

    bool connect_contour(
      std::vector<Segment_2>& contour) {

      bool success = false;

      success = clean_segments(contour);
      if (!success) return false;

      success = make_segments_collinear(contour);
      if (!success) return false;

      intersect_segments(contour);

      success = clean_and_intersect_segments(contour);
      if (success)
        success = clean_segments(contour);
      return success;
    }

    bool connect_polyline(
      std::vector<Segment_2>& contour) {

      bool success = false;

      /*
      success = clean_polyline(contour);
      if (!success) return false;

      success = make_segments_collinear(contour);
      if (!success) return false;

      intersect_polyline_segments(contour); */

      success = clean_and_intersect_polyline_segments(contour);
      return success;
    }

    bool clean_and_intersect_segments(
      std::vector<Segment_2>& contour) {

      const bool success = clean_segments(contour);
      if (!success) return false;

      intersect_segments(contour);

      for (const auto& segment : contour) {
        const auto& s = segment.source();
        const auto& t = segment.target();

        if (
          std::isnan(CGAL::to_double(s.x())) || std::isnan(CGAL::to_double(s.y())) ||
          std::isnan(CGAL::to_double(t.x())) || std::isnan(CGAL::to_double(t.y())) )
        return false;
      }
      return true;
    }

    bool clean_and_intersect_polyline_segments(
      std::vector<Segment_2>& contour) {

      const bool success = clean_polyline(contour);
      if (!success) return false;

      intersect_polyline_segments(contour);

      for (const auto& segment : contour) {
        const auto& s = segment.source();
        const auto& t = segment.target();

        if (
          std::isnan(CGAL::to_double(s.x())) || std::isnan(CGAL::to_double(s.y())) ||
          std::isnan(CGAL::to_double(t.x())) || std::isnan(CGAL::to_double(t.y())) )
        return false;
      }
      return true;
    }

    bool clean_segments(
      std::vector<Segment_2>& contour) {

      std::vector<Segment_2> clean;
      std::vector<Segment_2> segments;
      std::vector< std::vector<FT> > ratios;

      remove_zero_length_segments(contour, clean);
      if (clean.size() < 4)
        return false;

      const bool success = filter_out_wrong_segments(clean, segments, ratios);
      if (segments.size() < 4 || !success)
        return false;

      contour = segments;
      return true;
    }

    bool clean_polyline(
      std::vector<Segment_2>& contour) {

      std::vector<Segment_2> clean;
      std::vector<Segment_2> segments;
      std::vector< std::vector<FT> > ratios;

      remove_zero_length_segments(contour, clean);
      if (clean.size() < 1)
        return false;

      const bool success = filter_out_wrong_polyline_segments(
        clean, segments, ratios);
      if (segments.size() < 1 || !success)
        return false;

      contour = segments;
      return true;
    }

    bool optimize_contour(
      std::vector<Segment_2>& contour) {

      std::vector<Segment_2> clean;
      remove_zero_length_segments(contour, clean);
      if (clean.size() < 4)
        return false;

      std::map<std::size_t, std::size_t> seg_map;
      std::vector< std::vector<Segment_2> > groups;
      create_consecutive_groups(clean, groups, seg_map);

      /*
      std::size_t num_groups = 0;
      for (const auto& group : groups)
        if (group.size() > 1) ++num_groups;
      std::cout << "Num consecutive groups: " << num_groups << std::endl; */

      for (auto& group : groups)
        if (group.size() > 1)
          optimize_group(group);

      contour.clear();
      for (const auto& group : groups) {
        for (const auto& segment : group)
          contour.push_back(segment);
      }

      return true;
    }

    void create_consecutive_groups(
      const std::vector<Segment_2>& segments,
      std::vector< std::vector<Segment_2> >& groups,
      std::map<std::size_t, std::size_t>& seg_map) {

      groups.clear(); seg_map.clear();
      std::vector<bool> states(segments.size(), false);

      std::vector<Segment_2> group;
      std::size_t gr_idx = 0;

      const std::size_t num_segments = segments.size();
      for (std::size_t i = 0; i < num_segments; ++i) {
        const auto& segment_i = segments[i];
        if (states[i]) continue;

        group.clear();
        group.push_back(segment_i);
        seg_map[i] = gr_idx;
        states[i] = true;

        const std::size_t ip = (i + 1) % num_segments;
        if (ip != 0) {
          for (std::size_t j = ip; j < num_segments; ++j) {
            const auto& segment_j = segments[j];

            const FT angle   = angle_degree_2(segment_i, segment_j);
            const FT angle_2 = get_angle_2(angle);

            if (CGAL::abs(angle_2) <= m_angle_threshold) {
              group.push_back(segment_j); states[j] = true;
              seg_map[j] = gr_idx;
            } else break;
          }
        }
        groups.push_back(group);
        ++gr_idx;
      }
    }

    void optimize_group(
      std::vector<Segment_2>& segments) {

      std::vector<std::size_t> indices;
      indices.reserve(segments.size());
      for (std::size_t i = 0; i < segments.size(); ++i)
        indices.push_back(i);

      std::sort(indices.begin(), indices.end(),
      [&segments](const std::size_t idx1, const std::size_t idx2) -> bool {
        const FT length_1 = segments[idx1].squared_length();
        const FT length_2 = segments[idx2].squared_length();
        return length_1 > length_2;
      });

      std::vector< std::vector<Segment_2> > groups;
      std::map<std::size_t, std::size_t> seg_map;
      std::vector<bool> states(segments.size(), false);

      std::vector<Segment_2> group;
      std::size_t gr_idx = 0;

      for (std::size_t i = 0; i < indices.size(); ++i) {
        const int idxi = int(indices[i]);
        const auto& segment_i = segments[idxi];
        if (states[idxi]) continue;

        group.clear();
        group.push_back(segment_i);
        seg_map[idxi] = gr_idx;
        states[idxi] = true;

        const auto p =
        internal::middle_point_2(segment_i.source(), segment_i.target());

        const int idxip = idxi + 1;
        if (idxi < segments.size() - 1 && !states[idxip]) {
          for (int j = idxip; j < segments.size(); ++j) {
            const auto& segment_j = segments[j];
            if (states[j]) break;

            const Line_2 line = Line_2(segment_j.source(), segment_j.target());
            const auto q = line.projection(p);
            const FT distance = internal::distance(p, q);

            if (distance <= m_ordinate_bound) {
              group.push_back(segment_j); states[j] = true;
              seg_map[j] = gr_idx;
            } else break;
          }
        }

        const int idxim = idxi - 1;
        if (idxi > 0 && !states[idxim]) {
          auto j = idxim;
          while (j >= 0) {
            const auto& segment_j = segments[j];
            if (states[j]) break;

            const Line_2 line = Line_2(segment_j.source(), segment_j.target());
            const auto q = line.projection(p);
            const FT distance = internal::distance(p, q);

            if (distance <= m_ordinate_bound) {
              group.push_back(segment_j); states[j] = true;
              seg_map[j] = gr_idx;
            } else break;
            --j;
          }
        }

        groups.push_back(group);
        ++gr_idx;
      }

      std::vector<Line_2> lines;
      lines.reserve(groups.size());
      for (const auto& group : groups) {

        const Segment_2 segment = find_weighted_segment(group);
        const Line_2 line = Line_2(segment.source(), segment.target());
        lines.push_back(line);
      }

      std::vector<Segment_2> result;
      Point_2 p, q;
      for (std::size_t i = 0; i < segments.size(); ++i) {
        const std::size_t gr_idx = seg_map.at(i);
        const Line_2& line = lines[gr_idx];

        auto& segment = segments[i];

        const auto& s = segment.source();
        const auto& t = segment.target();

        p = line.projection(s);
        q = line.projection(t);
        segment = Segment_2(p, q);
      }

      for (std::size_t i = 0; i < segments.size(); ++i) {
        const std::size_t ip = (i + 1) % segments.size();
        if (ip == 0) {
          result.push_back(segments[i]);
          break;
        }

        const auto gri = seg_map.at(i);
        const auto grj = seg_map.at(ip);

        result.push_back(segments[i]);
        if (gri != grj) {

          const auto& segmenti = segments[i];
          const auto& segmentj = segments[ip];

          Line_2 line = Line_2(segmentj.source(), segmentj.target());
          auto source = internal::middle_point_2(segmenti.source(), segmenti.target());
          auto target = line.projection(source);
          Segment_2 orth = Segment_2(source, target);

          result.push_back(orth);
        }
      }
      segments = result;
    }

    void save_polylines(
      const std::vector<Segment_2>& segments,
      const std::string name) {

      CGAL_assertion(segments.size() > 0);
      std::vector< std::vector<Point_3> > polylines(segments.size());
      for (std::size_t i = 0; i < segments.size(); ++i) {
        const Point_2& s = segments[i].source();
        const Point_2& t = segments[i].target();

        polylines[i].push_back(Point_3(s.x(), s.y(), FT(0)));
        polylines[i].push_back(Point_3(t.x(), t.y(), FT(0)));
      }

      Saver<Traits> saver;
      saver.export_polylines(polylines, name);
    }

    bool make_segments_collinear(
      std::vector<Segment_2>& segments) {

      std::map<std::size_t, std::size_t> seg_map;
      std::vector< std::vector<Segment_2> > groups;
      create_collinear_groups(segments, groups, seg_map);

      std::size_t num_groups = 0;
      for (const auto& group : groups)
        if (group.size() > 1) ++num_groups;

      /* std::cout << "Num collinear groups: " << num_groups << std::endl; */

      std::vector<Line_2> lines;
      lines.reserve(groups.size());
      for (const auto& group : groups) {

        const Segment_2 segment = find_weighted_segment(group);
        const Line_2 line = Line_2(segment.source(), segment.target());
        lines.push_back(line);
      }

      for (std::size_t i = 0; i < segments.size(); ++i) {
        const std::size_t gr_idx = seg_map.at(i);
        const Line_2& line = lines[gr_idx];

        auto& segment = segments[i];

        const auto& s = segment.source();
        const auto& t = segment.target();

        // Do not use return here. It will not work! Tested!
        if (line.a() == FT(0) && line.b() == FT(0) && line.c() == FT(0))
          continue;

        const Point_2 p = line.projection(s);
        const Point_2 q = line.projection(t);

        segment = Segment_2(p, q);
      }
      return true;
    }

    void create_collinear_groups(
      const std::vector<Segment_2>& segments,
      std::vector< std::vector<Segment_2> >& groups,
      std::map<std::size_t, std::size_t>& seg_map) {

      groups.clear(); seg_map.clear();
      std::vector<bool> states(segments.size(), false);

      std::vector<Segment_2> group;
      std::size_t gr_idx = 0;

      for (std::size_t i = 0; i < segments.size(); ++i) {
        const auto& segment_i = segments[i];
        if (states[i]) continue;

        group.clear(); group.push_back(segment_i);
        seg_map[i] = gr_idx;
        states[i] = true;

        const auto p =
          internal::middle_point_2(segment_i.source(), segment_i.target());
        for (std::size_t j = 0; j < segments.size(); ++j) {
          const auto& segment_j = segments[j];
          if (states[j]) continue;

          const FT angle   = angle_degree_2(segment_i, segment_j);
          const FT angle_2 = get_angle_2(angle);

          if (CGAL::abs(angle_2) <= m_angle_threshold) {
            const Line_2 line = Line_2(segment_j.source(), segment_j.target());

            const auto q = line.projection(p);
            const FT distance = internal::distance(p, q);

            if (distance <= m_ordinate_bound) {
              group.push_back(segment_j); states[j] = true;
              seg_map[j] = gr_idx;
            }
          }
        }
        groups.push_back(group);
        ++gr_idx;
      }
    }

    void remove_zero_length_segments(
      const std::vector<Segment_2>& contour,
      std::vector<Segment_2>& segments) {

      segments.clear();
      for (const auto& segment : contour)
        if (segment.squared_length() > internal::tolerance<FT>())
          segments.push_back(segment);
    }

    bool filter_out_wrong_segments(
      const std::vector<Segment_2>& contour,
      std::vector<Segment_2>& segments,
      std::vector< std::vector<FT> >& ratios) {

      ratios.clear();
      segments.clear();

      const std::size_t n = contour.size();
      const std::size_t start = find_initial_index(contour);

      std::size_t i = start;
      std::vector<Segment_2> parallel_segments;
      std::size_t max_count = 0;
      do {

        const bool success = get_parallel_segments(
          contour, parallel_segments, i);
        if (!success) return false;

        Segment_2 segment;
        const FT sum_length =
        create_segment_from_parallel_segments(parallel_segments, segment);
        segments.push_back(segment);
        add_ratios(sum_length, parallel_segments, segment, ratios);

        ++max_count;
      } while (i != start && max_count < n * 2);
      if (max_count > n * 2) return false;
      return true;
    }

    bool filter_out_wrong_polyline_segments(
      const std::vector<Segment_2>& contour,
      std::vector<Segment_2>& segments,
      std::vector< std::vector<FT> >& ratios) {

      ratios.clear();
      segments.clear();

      const std::size_t n = contour.size();
      const std::size_t start = find_initial_polyline_index(contour);

      std::size_t i = start;
      std::vector<Segment_2> parallel_segments;
      std::size_t max_count = 0;
      do {

        const bool success = get_parallel_polyline_segments(
          contour, parallel_segments, i);
        if (!success) return false;

        Segment_2 segment;
        const FT sum_length =
        create_segment_from_parallel_segments(parallel_segments, segment);
        segments.push_back(segment);
        add_ratios(sum_length, parallel_segments, segment, ratios);

        ++max_count;
      } while (i != n && max_count < n);
      if (max_count > n) return false;
      return true;
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

        const auto pair = is_parallel_segment(sm, si, sp);
        const bool previous_is_orthogonal = !(pair.first);
        if (previous_is_orthogonal) return i;
      }
      return 0;
    }

    std::size_t find_initial_polyline_index(
      const std::vector<Segment_2>& contour) {

      /*
      const std::size_t n = contour.size();
      for (std::size_t i = 1; i < n - 1; ++i) {

        const std::size_t im = i - 1;
        const std::size_t ip = i + 1;

        const auto& si = contour[i];
        const auto& sm = contour[im];
        const auto& sp = contour[ip];

        const auto pair = is_parallel_segment(sm, si, sp);
        const bool previous_is_orthogonal = !(pair.first);
        if (previous_is_orthogonal) return i;
        if (!previous_is_orthogonal && i == 1) return 0;
      } */

      return 0;
    }

    std::pair<bool, bool> is_parallel_segment(
      const Segment_2& sm, const Segment_2& si, const Segment_2& sp) {

      const FT angle_mi = angle_degree_2(sm, si);
      const FT angle_mi_2 = get_angle_2(angle_mi);

      const FT angle_pi = angle_degree_2(si, sp);
      const FT angle_pi_2 = get_angle_2(angle_pi);

      const bool source_cond = ( CGAL::abs(angle_mi_2) <= m_angle_threshold );
      const bool target_cond = ( CGAL::abs(angle_pi_2) <= m_angle_threshold );

      return std::make_pair(source_cond, target_cond);
    }

    bool is_parallel_segment(
      const Segment_2& si, const Segment_2& sp) {

      const FT angle_pi = angle_degree_2(si, sp);
      const FT angle_pi_2 = get_angle_2(angle_pi);
      const bool target_cond = ( CGAL::abs(angle_pi_2) <= m_angle_threshold );
      return target_cond;
    }

    bool get_parallel_segments(
      const std::vector<Segment_2>& contour,
      std::vector<Segment_2>& parallel_segments,
      std::size_t& seed) {

      parallel_segments.clear();
      const std::size_t n = contour.size();

      std::size_t i = seed;
      bool next_is_parallel = false;
      std::size_t max_count = 0;
      do {

        const std::size_t im = (i + n - 1) % n;
        const std::size_t ip = (i + 1) % n;

        const auto& si = contour[i];
        const auto& sm = contour[im];
        const auto& sp = contour[ip];

        parallel_segments.push_back(si);
        const auto pair = is_parallel_segment(sm, si, sp);
        next_is_parallel = pair.second;
        i = ip;

        ++max_count;
      } while (next_is_parallel && max_count < n * 2);
      if (max_count > n * 2) return false;
      seed = i;
      return true;
    }

    bool get_parallel_polyline_segments(
      const std::vector<Segment_2>& contour,
      std::vector<Segment_2>& parallel_segments,
      std::size_t& seed) {

      parallel_segments.clear();
      const std::size_t n = contour.size();

      std::size_t i = seed;
      bool next_is_parallel = false;
      std::size_t max_count = 0;
      do {

        std::size_t im = std::size_t(-1);
        if (i > 0 && i <= n) im = i - 1;
        else im = std::size_t(-1);

        std::size_t ip = std::size_t(-1);
        if (i < n - 1 && i >= 0) ip = i + 1;
        else ip = std::size_t(-1);

        const auto& si = contour[i];
        parallel_segments.push_back(si);

        if (
          ip == std::size_t(-1)) {
          ++i; break;
        }

        if (
          im == std::size_t(-1)) {

          const auto& sp = contour[ip];
          next_is_parallel = is_parallel_segment(si, sp);
          i = ip;

        } else {

          const auto& sm = contour[im];
          const auto& sp = contour[ip];
          const auto pair = is_parallel_segment(sm, si, sp);
          next_is_parallel = pair.second;
          i = ip;
        }

        ++max_count;
      } while (next_is_parallel && max_count < n);
      if (max_count > n) return false;
      seed = i;
      return true;
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

      source = Point_2(x1, y1);
      target = Point_2(x2, y2);

      if (source == target)
        return find_longest_segment(segments);
      return Segment_2(source, target);
    }

    FT create_segment_from_parallel_segments(
      const std::vector<Segment_2>& parallel_segments,
      Segment_2& result) {

      Segment_2 ref_segment = find_weighted_segment(parallel_segments);
      const Line_2 line =
      Line_2(ref_segment.source(), ref_segment.target());

      FT sum_length = FT(0);
      std::vector<Point_2> points;
      for (const auto& segment : parallel_segments) {

        const Point_2 p = line.projection(segment.source());
        const Point_2 q = line.projection(segment.target());

        points.push_back(p);
        points.push_back(q);

        sum_length += internal::distance(p, q);
      }
      update_segment(points, ref_segment);
      result = ref_segment;
      return sum_length;
    }

    Segment_2 find_weighted_segment(
      const std::vector<Segment_2>& segments) {

      std::vector<FT> weights;
      compute_distance_weights(segments, weights);
      const Segment_2 ref_segment = find_central_segment(segments);
      const Segment_2 result =
        compute_weighted_segment(segments, weights, ref_segment);

      if (result.source() == result.target())
        return ref_segment;
      return result;
    }

    void compute_distance_weights(
      const std::vector<Segment_2>& segments,
      std::vector<FT>& weights) {

      weights.clear();
      weights.reserve(segments.size());

      FT sum_distance = FT(0);
      for (const auto& segment : segments) {
        const FT distance =
          internal::distance(segment.source(), segment.target());
        sum_distance += distance;

        weights.push_back(distance);
      }

      for (auto& weight : weights)
        weight /= sum_distance;
    }

    Segment_2 compute_weighted_segment(
      const std::vector<Segment_2>& segments,
      const std::vector<FT>& weights,
      const Segment_2& ref_segment) {

      const Point_2& s = ref_segment.source();
      const Point_2& t = ref_segment.target();

      const Point_2 b = internal::middle_point_2(s, t);

      Vector_2 dir = Vector_2(FT(0), FT(0));
      for (std::size_t i = 0; i < weights.size(); ++i) {
        const FT weight = weights[i];

        const Segment_2& segment = segments[i];
        const Line_2 line = Line_2(segment.source(), segment.target());
        const Point_2 p = line.projection(b);

        const Vector_2 v = Vector_2(b, p);
        dir += v * weight;
      }

      const Point_2 news = s + dir;
      const Point_2 newt = t + dir;

      return Segment_2(news, newt);
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

    void add_ratios(
      const FT sum_length,
      const std::vector<Segment_2>& parallel_segments,
      const Segment_2& segment,
      std::vector< std::vector<FT> >& ratios) {

      const FT ref_length =
      internal::distance(segment.source(), segment.target());
      const FT error = (FT(1) - sum_length / ref_length) /
        static_cast<FT>(parallel_segments.size());

      FT length = FT(0);
      std::vector<FT> ds;
      ds.push_back(length);

      for (std::size_t k = 0; k < parallel_segments.size(); ++k) {
        length += internal::distance(
          parallel_segments[k].source(), parallel_segments[k].target());
        const FT d = length / ref_length + static_cast<FT>(k + 1) * error;
        ds.push_back(d);
      }
      ratios.push_back(ds);
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

    void intersect_polyline_segments(
      std::vector<Segment_2>& segments) {

      const std::size_t n = segments.size();
      for (std::size_t i = 1; i < n - 1; ++i) {

        const std::size_t im = i - 1;
        const std::size_t ip = i + 1;

        auto& si = segments[i];
        const auto& sm = segments[im];
        const auto& sp = segments[ip];

        intersect_segment(sm, si, sp);
      }
    }

    void intersect_segment(
      const Segment_2& sm, Segment_2& si, const Segment_2& sp) {

      Point_2 source = si.source();
      Point_2 target = si.target();

      const Line_2 line_1 = Line_2(sm.source(), sm.target());
      const Line_2 line_2 = Line_2(si.source(), si.target());
      const Line_2 line_3 = Line_2(sp.source(), sp.target());

      const bool success1 = intersect_2(line_1, line_2, source);
      const bool success2 = intersect_2(line_2, line_3, target);

      if (!success1) source = si.source();
      if (!success2) target = si.target();

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

    void split_segments(
      const std::vector<Segment_2>& segments,
      const std::vector< std::vector<FT> >& ratios,
      std::vector<Segment_2>& splitted) {

      splitted.clear();
      for (std::size_t i = 0; i < ratios.size(); ++i) {

        const Segment_2& segment = segments[i];
        const FT ref_length = internal::distance(segment.source(), segment.target());
        Vector_2 direction = Vector_2(segment.source(), segment.target());
        internal::normalize(direction);
        const Vector_2 start = Vector_2(
          Point_2(FT(0), FT(0)), segment.source());

        for (std::size_t j = 0; j < ratios[i].size() - 1; ++j) {
          const std::size_t jp = j + 1;

          const Vector_2 end1 = start + direction * (ratios[i][j]  * ref_length);
          const Vector_2 end2 = start + direction * (ratios[i][jp] * ref_length);

          const Point_2 source = Point_2(end1.x(), end1.y());
          const Point_2 target = Point_2(end2.x(), end2.y());

          const Segment_2 seg = Segment_2(source, target);
          splitted.push_back(seg);
        }
      }
    }
  };

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_INTERNAL_POLYGON_REGULARIZER_H
