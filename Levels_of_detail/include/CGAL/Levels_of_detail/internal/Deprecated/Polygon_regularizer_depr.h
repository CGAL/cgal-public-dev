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

#ifndef CGAL_LEVELS_OF_DETAIL_INTERNAL_POLYGON_REGULARIZER_DEPR_H
#define CGAL_LEVELS_OF_DETAIL_INTERNAL_POLYGON_REGULARIZER_DEPR_H

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
  class Polygon_regularizer_depr {

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

    Polygon_regularizer_depr(
      const FT min_length,
      const FT angle_bound) :
    m_min_length(min_length),
    m_angle_bound(angle_bound),
    m_pi(static_cast<FT>(CGAL_PI)),
    m_angle_threshold(FT(5)),
    m_bound_min(m_angle_bound / FT(3)),
    m_bound_max(FT(90) - m_bound_min)
    { }

    void compute_longest_direction(
      const std::vector< std::vector<Segment_2> >& contours) {

      m_bounds.clear();
      m_bounds.resize(1);
      m_bounds[0] = std::make_pair(FT(45), FT(45));

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

    void compute_longest_directions(
      const std::vector< std::vector<Segment_2> >& contours) {

      m_bounds.clear();
      m_bounds.resize(contours.size());

      for (std::size_t k = 0; k < contours.size(); ++k)
        m_bounds[k] = std::make_pair(FT(45), FT(45));

      m_skip.clear();
      m_skip.resize(contours.size());
      for (std::size_t k = 0; k < contours.size(); ++k)
        m_skip[k] = std::make_pair(k, find_longest_segment(contours[k]));

      m_longest.clear();
      m_longest.resize(contours.size());
      for (std::size_t k = 0; k < contours.size(); ++k)
        m_longest[k] = contours[m_skip[k].first][m_skip[k].second];

      m_groups.clear();
      m_groups.resize(contours.size());
      for (std::size_t k = 0; k < contours.size(); ++k) {
        m_groups[k].resize(contours[k].size());
        for (std::size_t i = 0; i < contours[k].size(); ++i) {
          m_groups[k][i] = k;
        }
      }
    }

    void compute_multiple_directions(
      const std::vector< std::vector<Segment_2> >& cnt) {

      std::vector< std::pair<Segment_2, bool> > segments;
      std::vector< std::vector< std::pair<Segment_2, bool> > > contours;

      for (std::size_t k = 0; k < cnt.size(); ++k) {
        const auto& contour = cnt[k];
        segments.clear();
        for (std::size_t i = 0; i < contour.size(); ++i) {
          if (internal::distance(contour[i].source(), contour[i].target()) > m_min_length * FT(2))
            segments.push_back(std::make_pair(contour[i], true));
          else
            segments.push_back(std::make_pair(contour[i], false));
        }
        contours.push_back(segments);
      }

      std::vector<Size_pair> input;
      for (std::size_t k = 0; k < contours.size(); ++k)
        for (std::size_t i = 0; i < contours[k].size(); ++i)
          input.push_back(std::make_pair(k, i));

      sort_input(contours, input);
      std::vector<bool> states(input.size(), false);

      m_bounds.clear(); m_skip.clear();
      m_longest.clear(); m_groups.clear();

      m_groups.resize(contours.size());
      for (std::size_t k = 0; k < contours.size(); ++k) {
        m_groups[k].resize(contours[k].size());
        for (std::size_t i = 0; i < contours[k].size(); ++i) {
          m_groups[k][i] = std::size_t(-1);
        }
      }

      bool apply = true; std::size_t gr_idx = 0;
      do {
        apply = get_next_direction(
          contours, input, gr_idx, states);
        ++gr_idx;
      } while (apply);

      if (m_longest.size() == 0) {
        compute_longest_direction(cnt);

        std::cout << "Num outer directions: " << m_longest.size() << std::endl;
        return;
      }

      std::cout << "Num outer directions: " << m_longest.size() << std::endl;
      for (std::size_t k = 0; k < contours.size(); ++k) {
        for (std::size_t i = 0; i < contours[k].size(); ++i) {
          if (m_groups[k][i] == std::size_t(-1)) {

            const std::size_t m = contours[k].size();
            std::size_t im = (i + m - 1) % m;
            std::size_t ip = (i + 1) % m;

            bool stop = false;
            std::size_t max_count = 0;
            do {
              if (contours[k][im].second) {
                m_groups[k][i] = m_groups[k][im];
                break;
              }
              if (contours[k][ip].second) {
                m_groups[k][i] = m_groups[k][ip];
                break;
              }
              im = (im + m - 1) % m;
              ip = (ip + 1) % m;

              if (im == i || ip == i)
                stop = true;

              ++max_count;
            } while (!stop && max_count < m * 2);

            if (stop || max_count >= m * 2) {
              m_groups[k][i] = 0;
            }
          }
        }
      }
    }

    bool get_next_direction(
      const std::vector< std::vector< std::pair<Segment_2, bool> > >& contours,
      const std::vector<Size_pair>& input,
      const std::size_t gr_idx,
      std::vector<bool>& states) {

      // Add new group.
      std::size_t longest_idx = std::size_t(-1);
      for (std::size_t i = 0; i < states.size(); ++i) {
        if (!states[i] && contours[input[i].first][input[i].second].second) {
          longest_idx = i;
          break;
        }
      }
      if (longest_idx == std::size_t(-1))
        return false;

      const auto& longest_pair = input[longest_idx];
      const Segment_2& longest =
        contours[longest_pair.first][longest_pair.second].first;

      // Fill in groups.
      for (std::size_t i = 0; i < states.size(); ++i) {
        if (i == longest_idx) {
          m_groups[longest_pair.first][longest_pair.second] = gr_idx;
          states[i] = true;
          continue;
        }

        if (!states[i]) {
          if (contours[input[i].first][input[i].second].second) {
            const auto& pair = input[i];
            const auto& si = contours[pair.first][pair.second].first;

            const FT angle = angle_degree_2(longest, si);
            const FT angle_2 = get_angle_2(angle);

            if (
              (CGAL::abs(angle_2) <= m_bound_min) ||
              (CGAL::abs(angle_2) >= m_bound_max) )  {

              m_groups[pair.first][pair.second] = gr_idx;
              states[i] = true;
              continue;
            }
          }
        }
      }

      m_bounds.push_back(std::make_pair(FT(45), FT(45)));
      m_skip.push_back(longest_pair);
      m_longest.push_back(longest);

      return true;
    }

    void sort_input(
      const std::vector< std::vector< std::pair<Segment_2, bool> > >& contours,
      std::vector<Size_pair>& input) {

      std::sort(input.begin(), input.end(),
      [&contours](const Size_pair& a, const Size_pair& b) -> bool {
        const FT length_1 = (contours[a.first][a.second]).first.squared_length();
        const FT length_2 = (contours[b.first][b.second]).first.squared_length();
        return length_1 > length_2;
      });
    }

    void regularize_contours(
      std::vector< std::vector<Segment_2> >& contours) {

      if (m_angle_bound == FT(0))
        return;

      std::vector< std::vector<Segment_2> > initials, finals;
      /* std::vector< std::vector< std::pair<std::vector<Point_2>, FT> > > points; */

      /* std::vector<Segment_2> clean; */
      initials.reserve(contours.size());
      for (std::size_t k = 0; k < contours.size(); ++k) {
        const auto& contour = contours[k];

        /*
        remove_zero_length_segments(contour, clean);
        if (clean.size() < 4)
          continue;
        initials.push_back(clean); */

        initials.push_back(contour);

        /* points.push_back(contour_points[k]); */
      }

      for (std::size_t k = 0; k < initials.size(); ++k) {
        auto contour = initials[k];

        rotate_contour(k, contour);
        correct_contour(k, contour);
        const bool success = connect_contour(contour, false);
        if (success)
          finals.push_back(contour);
      }

      contours = finals;

      /*
      std::vector< std::vector<Segment_2> > result;
      apply_least_squares(initials, finals, points, result);

      for (std::size_t k = 0; k < result.size(); ++k) {
        auto contour = result[k];
        connect_contour(contour, false);
        contours.push_back(contour);
      } */
    }

    void apply_least_squares(
      const std::vector< std::vector<Segment_2> >& initials,
      const std::vector< std::vector<Segment_2> >& finals,
      const std::vector< std::vector< std::pair<std::vector<Point_2>, FT> > >& points,
      std::vector< std::vector<Segment_2> >& result) {

      result.clear();
      result.resize(initials.size());
      for (std::size_t k = 0; k < initials.size(); ++k) {
        const auto& inp_segments = initials[k];
        const auto& fin_segments = finals[k];
        const auto& pts = points[k];

        for (std::size_t kk = 0; kk < inp_segments.size(); ++kk) {
          const FT error = compute_error(fin_segments[kk], pts[kk].first);
          if ( ( error <= pts[kk].second ) ||
               ( internal::distance(inp_segments[kk].source(), inp_segments[kk].target()) <= m_min_length ) ) {
            result[k].push_back(fin_segments[kk]);
          } else {
            result[k].push_back(inp_segments[kk]);
          }
        }
      }
    }

    FT compute_error(
      const Segment_2& segment,
      const std::vector<Point_2>& points) {

      if (points.size() == 0)
        return FT(0);

      const Line_2 line = Line_2(segment.source(), segment.target());
      FT dist = FT(0);
      for (const auto& p : points) {
        const auto& q = line.projection(p);
        dist += internal::distance(p, q);
      }
      CGAL_assertion(points.size() > 0);
      dist /= static_cast<FT>(points.size());
      return dist;
    }

    void rotate_contour(
      const std::size_t k,
      std::vector<Segment_2>& contour) {

      for (std::size_t i = 0; i < contour.size(); ++i) {
        const std::size_t gr_idx = m_groups[k][i];
        if (k == m_skip[gr_idx].first && i == m_skip[gr_idx].second)
          continue;
        if (gr_idx == std::size_t(-1))
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
        if (gr_idx == std::size_t(-1))
          continue;

        const std::size_t im = (i + n - 1) % n;
        const std::size_t ip = (i + 1) % n;

        auto& si = contour[i];
        const auto& sm = contour[im];
        const auto& sp = contour[ip];

        const FT length = internal::distance(si.source(), si.target());
        if (length <= m_min_length)
          correct_segment(sm, si, sp);
      }
    }

    bool connect_contour(
      std::vector<Segment_2>& contour,
      const bool split) {

      std::vector<Segment_2> clean;
      std::vector<Segment_2> segments;
      std::vector<Segment_2> splitted;

      std::vector< std::vector<FT> > ratios;

      remove_zero_length_segments(contour, clean);
      if (clean.size() < 4)
        return false;

      const bool success = filter_out_wrong_segments(clean, segments, ratios);
      if (segments.size() < 4 || !success)
        return false;

      intersect_segments(segments);

      if (split) {
        split_segments(segments, ratios, splitted);
        contour = splitted;
        return true;
      }
      contour = segments;
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

    bool filter_out_wrong_segments(
      const std::vector<Segment_2>& contour,
      std::vector<Segment_2>& segments,
      std::vector< std::vector<FT> >& ratios) {

      ratios.clear();
      segments.clear();
      const std::size_t n = contour.size();
      const std::size_t start = find_initial_index(contour);
      std::size_t i = start;
      std::vector< std::vector<Segment_2> > tmp(1);
      std::size_t max_count_out = 0;
      do {

        tmp[0].clear(); std::size_t j = i;
        bool next_is_parallel = false;
        std::size_t max_count_in = 0;
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

          ++max_count_in;
        } while (next_is_parallel && max_count_in < n * 2);

        if (max_count_in >= n * 2)
          return false;

        /*
        const auto data = find_longest_segment(tmp);
        segments.push_back(tmp[data.first][data.second]); */

        auto segment = find_central_segment(tmp[0]);
        const Line_2 line = Line_2(segment.source(), segment.target());

        FT sum_length = FT(0);
        std::vector<Point_2> points;
        for (auto& seg : tmp[0]) {
          const Point_2 p = line.projection(seg.source());
          const Point_2 q = line.projection(seg.target());
          points.push_back(p);
          points.push_back(q);
          seg = Segment_2(p, q);
          sum_length += internal::distance(p, q);
        }
        update_segment(points, segment);
        segments.push_back(segment);
        const FT ref_length =
        internal::distance(segment.source(), segment.target());

        FT length = FT(0);
        std::vector<FT> ds;
        ds.push_back(length);

        const FT error = (FT(1) - sum_length / ref_length) /
          static_cast<FT>(tmp[0].size());

        for (std::size_t kk = 0; kk < tmp[0].size(); ++kk) {
          length += internal::distance(tmp[0][kk].source(), tmp[0][kk].target());
          const FT d = length / ref_length + (kk + 1) * error;
          ds.push_back(d);
        }
        ratios.push_back(ds);
        i = j;

        ++max_count_out;
      } while (i != start && max_count_out < n * 2);

      if (max_count_out >= n * 2)
        return false;
      return true;
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
      return 0;
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
    const FT m_bound_min, m_bound_max;

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

    std::size_t find_longest_segment(
      const std::vector<Segment_2>& contour) {

      std::size_t seg_idx = std::size_t(-1);
      FT max_length = -FT(1);
      for (std::size_t i = 0; i < contour.size(); ++i) {

        const FT length = contour[i].squared_length();
        if (length > max_length) {

          max_length = length;
          seg_idx = i;
        }
      }
      return seg_idx;
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

      const bool success1 = intersect_2(line_1, line_2, source);
      const bool success2 = intersect_2(line_2, line_3, target);

      if (!success1) source = si.source();
      if (!success2) target = si.target();

      /*
      source = line_1.projection(si.source());
      target = line_3.projection(si.target()); */

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

#endif // CGAL_LEVELS_OF_DETAIL_INTERNAL_POLYGON_REGULARIZER_DEPR_H
