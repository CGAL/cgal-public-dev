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
#include <set>
#include <vector>
#include <utility>
#include <algorithm>

// CGAL includes.
#include <CGAL/assertions.h>
#include <CGAL/property_map.h>
#include <CGAL/point_generators_2.h>

// Internal includes.
#include <CGAL/Levels_of_detail/internal/utils.h>
#include <CGAL/Levels_of_detail/internal/struct.h>

// Spatial search.
#include <CGAL/Levels_of_detail/internal/Spatial_search/K_neighbor_query.h>

// Testing.
#include "../../../../../test/Levels_of_detail/include/Saver.h"

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

  template<typename GeomTraits>
  class Segment_regularizer {

  public:
    using Traits = GeomTraits;

    using FT = typename Traits::FT;
    using Point_2 = typename Traits::Point_2;
    using Point_3 = typename Traits::Point_3;
    using Segment_2 = typename Traits::Segment_2;
    using Vector_2 = typename Traits::Vector_2;
    using Line_2 = typename Traits::Line_2;

    using FT_pair = std::pair<FT, FT>;
    using Size_pair = std::pair<std::size_t, std::size_t>;
    using Indices = std::vector<std::size_t>;
    using Seg_pair = std::pair<Segment_2, bool>;

    struct Range_data {
      std::size_t seg_i = std::size_t(-1);
      std::size_t seg_j = std::size_t(-1);

      bool is_corner = false;
      std::size_t gr_idx = std::size_t(-1);
    };

    using Point_pair = std::pair<Point_2, Range_data>;
    using Point_map = CGAL::First_of_pair_property_map<Point_pair>;
    using K_neighbor_query =
    internal::K_neighbor_query<Traits, std::vector<Point_pair>, Point_map>;
    using Point_generator = CGAL::Points_on_segment_2<Point_2>;

    Segment_regularizer(
      const FT min_length,
      const FT angle_bound,
      const FT angle_threshold = FT(5)) :
    m_min_length(min_length),
    m_angle_bound(angle_bound),
    m_angle_threshold(angle_threshold),
    m_bound_min(m_angle_bound),
    m_bound_max(FT(90) - m_bound_min),
    m_pi(static_cast<FT>(CGAL_PI)),
    m_k(10),
    m_num_samples_per_segment(m_k * 2)
    { }

    void compute_longest_direction(
      const std::vector<Segment_2>& outer_segments,
      const std::vector< std::vector<Segment_2> >& contours) {

      m_bounds.clear();
      m_bounds.resize(1);
      m_bounds[0] = std::make_pair(FT(45), FT(45));

      const std::size_t seg_idx = find_longest_segment(outer_segments);

      m_longest.clear();
      m_longest.resize(1);
      m_longest[0] = outer_segments[seg_idx];

      make_default_groups(contours, 0, m_groups);
    }

    void compute_multiple_directions(
      const std::vector<Segment_2>& segments_outer,
      const std::vector< std::vector<Segment_2> >& contours) {

      std::vector< std::vector<Seg_pair> > contours_outer;
      create_contours_from_segments(segments_outer, contours_outer);

      std::vector<FT_pair> bounds_outer;
      std::vector<Size_pair> skip_outer;
      std::vector<Segment_2> longest_outer;
      std::vector<Indices> groups_outer;

      get_multiple_directions(
        contours_outer,
        bounds_outer, skip_outer, longest_outer, groups_outer);

      unify_along_contours(segments_outer, contours_outer, groups_outer);
      std::cout << "Num outer directions: " << longest_outer.size() << std::endl;

      std::vector<Point_pair> pair_range;
      create_pair_range(contours_outer, groups_outer, pair_range);

      m_bounds.clear();
      m_bounds = bounds_outer;

      m_longest.clear();
      m_longest = longest_outer;

      m_groups.clear();
      make_default_groups(contours, std::size_t(-1), m_groups);
      assign_groups_using_kd_tree(
        contours, pair_range, m_longest, m_groups);
    }

    void regularize_contours(
      std::vector< std::vector<Segment_2> >& contours) {

      if (m_angle_bound == FT(0))
        return;

      m_saved = contours;
      for (std::size_t k = 0; k < contours.size(); ++k) {
        auto& contour = contours[k];
        rotate_contour(k, contour);
        if (contour.size() >= 3)
          correct_contour(k, contour);
      }
    }

  private:
    const FT m_min_length;
    const FT m_angle_bound;
    const FT m_angle_threshold;

    const FT m_bound_min, m_bound_max;
    const FT m_pi;

    const std::size_t m_k;
    const std::size_t m_num_samples_per_segment;

    std::vector<FT_pair> m_bounds;
    std::vector<Segment_2> m_longest;
    std::vector<Indices> m_groups;

    std::vector< std::vector<Segment_2> > m_saved;

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

    template<typename T>
    void make_default_groups(
      const std::vector< std::vector<T> >& contours,
      const std::size_t value,
      std::vector<Indices>& groups) {

      groups.clear(); groups.resize(contours.size());
      for (std::size_t k = 0; k < contours.size(); ++k)
        groups[k].resize(contours[k].size(), value);
    }

    void create_contours_from_segments(
      const std::vector<Segment_2>& segments,
      std::vector< std::vector<Seg_pair> >& contours) {

      contours.clear();
      std::vector<Seg_pair> contour;
      for (const auto& segment : segments) {
        const auto& s = segment.source();
        const auto& t = segment.target();

        if (internal::distance(s, t) >= m_min_length * FT(2))
          contour.push_back(std::make_pair(segment, true));
        else
          contour.push_back(std::make_pair(segment, false));
      }
      contours.push_back(contour);
    }

    void get_multiple_directions(
      const std::vector< std::vector<Seg_pair> >& contours,
      std::vector<FT_pair>& bounds,
      std::vector<Size_pair>& skip,
      std::vector<Segment_2>& longest,
      std::vector<Indices>& groups) {

      make_default_groups(contours, std::size_t(-1), groups);

      std::vector<Size_pair> input;
      for (std::size_t k = 0; k < contours.size(); ++k)
        for (std::size_t i = 0; i < contours[k].size(); ++i)
          input.push_back(std::make_pair(k, i));

      sort_input(contours, input);
      std::vector<bool> states(input.size(), false);

      bool apply = true; std::size_t gr_idx = 0;
      do {
        apply = get_next_direction(
          contours, input, gr_idx, states,
          bounds, skip, longest, groups);
          ++gr_idx;
      } while (apply);

      if (longest.size() == 0) {

        bounds.push_back(std::make_pair(FT(45), FT(45)));
        const auto longest_pair = find_longest_segment(contours);
        skip.push_back(longest_pair);
        const auto& longest_segment =
          (contours[skip[0].first][skip[0].second]).first;
        longest.push_back(longest_segment);
        make_default_groups(contours, 0, groups);
      }
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
      std::vector<bool>& states,
      std::vector<FT_pair>& bounds,
      std::vector<Size_pair>& skip,
      std::vector<Segment_2>& longest,
      std::vector<Indices>& groups) {

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
          groups[longest_pair.first][longest_pair.second] = gr_idx;
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

            groups[pair.first][pair.second] = gr_idx;
            states[i] = true; continue;
          }
        }
      }

      longest.push_back(longest_segment);
      bounds.push_back(std::make_pair(FT(45), FT(45)));
      skip.push_back(longest_pair);

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

    void unify_along_contours(
      const std::vector<Segment_2>& segments,
      const std::vector< std::vector<Seg_pair> >& contours,
      std::vector<Indices>& groups) {

      Range_data data;
      for (std::size_t k = 0; k < contours.size(); ++k) {
        for (std::size_t i = 0; i < contours[k].size(); ++i) {
          if (contours[k][i].second)
            continue;

          if (groups[k][i] == std::size_t(-1)) {
            const auto& segment = contours[k][i].first;
            find_corner(segment, segments, data);

            const std::size_t m = contours[k].size();

            std::size_t im = data.seg_i;
            std::size_t ip = data.seg_j;

            bool stop = false;
            std::size_t max_count = 0;
            do {

              if (im != std::size_t(-1) && contours[k][im].second) {
                groups[k][i] = groups[k][im]; break;
              }

              if (ip != std::size_t(-1) && contours[k][ip].second) {
                groups[k][i] = groups[k][ip]; break;
              }

              if (im != std::size_t(-1)) {
                const auto& sm = contours[k][im].first;
                find_corner(sm, segments, data);
                im = data.seg_i;
              }

              if (ip != std::size_t(-1)) {
                const auto& sp = contours[k][ip].first;
                find_corner(sp, segments, data);
                ip = data.seg_j;
              }

              if (im == i || ip == i) stop = true;
              ++max_count;

            } while (!stop && max_count < m * 2);
            if (stop || max_count >= m * 2)
              groups[k][i] = 0;
          }
        }
      }
    }

    void find_corner(
      const Segment_2& ref,
      const std::vector<Segment_2>& segments_outer,
      Range_data& data) {

      const auto& p = ref.source();
      const auto& q = ref.target();

      data.seg_i = std::size_t(-1);
      data.seg_j = std::size_t(-1);

      bool found_1 = false, found_2 = false;
      const std::size_t n = segments_outer.size();
      for (std::size_t i = 0; i < n; ++i) {
        const auto& segment = segments_outer[i];

        const auto& s = segment.source();
        const auto& t = segment.target();

        if (p == t) {
          found_1 = true; data.seg_i = i;
          if (found_2) break;
          else continue;
        }

        if (q == s) {
          found_2 = true; data.seg_j = i;
          if (found_1) break;
          else continue;
        }
      }
    }

    void create_pair_range(
      const std::vector< std::vector<Seg_pair> >& contours,
      const std::vector<Indices>& groups,
      std::vector<Point_pair>& pair_range) {

      pair_range.clear();
      std::vector<Point_2> samples;
      for (std::size_t k = 0; k < contours.size(); ++k) {
        for (std::size_t i = 0; i < contours[k].size(); ++i) {
          const auto& segment = contours[k][i].first;

          const auto& s = segment.source();
          const auto& t = segment.target();

          samples.clear();
          Point_generator generator(s, t, m_num_samples_per_segment);
          std::copy_n(generator, m_num_samples_per_segment - 1,
          std::back_inserter(samples));

          const std::size_t gr_idx = groups[k][i];
          for (const auto& p : samples) {
            Range_data data;
            data.gr_idx = gr_idx;
            pair_range.push_back(std::make_pair(p, data));
          }
        }
      }
    }

    void assign_groups_using_kd_tree(
      const std::vector< std::vector<Segment_2> >& contours,
      const std::vector<Point_pair>& pair_range,
      const std::vector<Segment_2>& longest,
      std::vector<Indices>& groups) {

      Point_map point_map;
      K_neighbor_query neighbor_query(pair_range, FT(m_k), point_map);

      Indices neighbors;
      for (std::size_t k = 0; k < contours.size(); ++k) {
        for (std::size_t l = 0; l < contours[k].size(); ++l) {
          const auto& segment = contours[k][l];

          std::size_t gr_idx_i = std::size_t(-1);
          std::size_t gr_idx_j = std::size_t(-1);
          FT abs_angle_i_2 = FT(0), abs_angle_j_2 = FT(0);

          get_best_angles(segment, pair_range, longest,
          neighbor_query, neighbors,
          gr_idx_i, gr_idx_j, abs_angle_i_2, abs_angle_j_2);

          if (
            abs_angle_i_2 <= FT(45) && abs_angle_j_2 <= FT(45) &&
            abs_angle_i_2 <= abs_angle_j_2) {
            groups[k][l] = gr_idx_i;
            continue;
          }

          if (
            abs_angle_i_2 <= FT(45) && abs_angle_j_2 <= FT(45) &&
            abs_angle_i_2 > abs_angle_j_2) {
            groups[k][l] = gr_idx_j;
            continue;
          }

          if (
            abs_angle_i_2 > FT(45) && abs_angle_j_2 > FT(45) &&
            abs_angle_i_2 > abs_angle_j_2) {
            groups[k][l] = gr_idx_i;
            continue;
          }

          if (
            abs_angle_i_2 > FT(45) && abs_angle_j_2 > FT(45) &&
            abs_angle_i_2 <= abs_angle_j_2) {
            groups[k][l] = gr_idx_j;
            continue;
          }

          if (
            abs_angle_i_2 <= FT(45) && abs_angle_j_2 > FT(45) &&
            abs_angle_i_2 < ( FT(90) - abs_angle_j_2 ) ) {
            groups[k][l] = gr_idx_i;
            continue;
          }

          if (
            abs_angle_i_2 <= FT(45) && abs_angle_j_2 > FT(45) &&
            abs_angle_i_2 >= ( FT(90) - abs_angle_j_2 ) ) {
            groups[k][l] = gr_idx_j;
            continue;
          }

          if (
            abs_angle_i_2 > FT(45) && abs_angle_j_2 <= FT(45) &&
            abs_angle_j_2 < ( FT(90) - abs_angle_i_2 ) ) {
            groups[k][l] = gr_idx_j;
            continue;
          }

          if (
            abs_angle_i_2 > FT(45) && abs_angle_j_2 <= FT(45) &&
            abs_angle_j_2 >= ( FT(90) - abs_angle_i_2 ) ) {
            groups[k][l] = gr_idx_i;
            continue;
          }
          groups[k][l] = 0;
        }
      }
    }

    void get_best_angles(
      const Segment_2& segment,
      const std::vector<Point_pair>& pair_range,
      const std::vector<Segment_2>& longest,
      K_neighbor_query& neighbor_query,
      Indices& neighbors,
      std::size_t& gr_idx_i, std::size_t& gr_idx_j,
      FT& abs_angle_i_2, FT& abs_angle_j_2) {

      const auto& s = segment.source();
      const auto& t = segment.target();

      gr_idx_i = get_group_index(
        s, segment, pair_range, longest, neighbor_query, neighbors);
      gr_idx_j = get_group_index(
        t, segment, pair_range, longest, neighbor_query, neighbors);

      if (gr_idx_i == std::size_t(-1) || gr_idx_j == std::size_t(-1))
        return;

      const auto& longest_i = longest[gr_idx_i];
      const auto& longest_j = longest[gr_idx_j];

      const FT angle_i   = angle_degree_2(segment, longest_i);
      const FT angle_i_2 = get_angle_2(angle_i);

      const FT angle_j   = angle_degree_2(longest_j, segment);
      const FT angle_j_2 = get_angle_2(angle_j);

      abs_angle_i_2 = CGAL::abs(angle_i_2);
      abs_angle_j_2 = CGAL::abs(angle_j_2);
    }

    std::size_t get_group_index(
      const Point_2& query,
      const Segment_2& segment,
      const std::vector<Point_pair>& pair_range,
      const std::vector<Segment_2>& longest,
      K_neighbor_query& neighbor_query,
      Indices& neighbors) {

      neighbor_query(query, neighbors);
      std::size_t idx = std::size_t(-1);
      FT angle_min = internal::max_value<FT>();

      for (const std::size_t neighbor : neighbors) {
        const std::size_t gr_idx = pair_range[neighbor].second.gr_idx;

        const FT angle   = angle_degree_2(longest[gr_idx], segment);
        const FT angle_2 = get_angle_2(angle);

        const FT abs_angle_2 = CGAL::abs(angle_2);
        if (abs_angle_2 <= FT(45) && abs_angle_2 < angle_min ) {
          angle_min = abs_angle_2;
          idx = neighbor;
        }

        if (abs_angle_2 >= FT(45) && angle_min > ( FT(90) - abs_angle_2 ) ) {
          angle_min = FT(90) - abs_angle_2;
          idx = neighbor;
        }
      }
      if (idx == std::size_t(-1)) return idx;
      return pair_range[idx].second.gr_idx;
    }

    void rotate_contour(
      const std::size_t k,
      std::vector<Segment_2>& contour) {

      for (std::size_t i = 0; i < contour.size(); ++i) {
        const std::size_t gr_idx = m_groups[k][i];
        if (gr_idx == std::size_t(-1))
          continue;

        auto& segment = contour[i];
        const auto& longest_segment = m_longest[gr_idx];
        const auto& bounds = m_bounds[gr_idx];

        const bool success = rotate_segment(longest_segment, bounds, segment);
        if (!success)
          m_groups[k][i] = std::size_t(-1);
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
      segments.push_back(contour[0]);
      for (std::size_t i = 1; i < n - 1; ++i) {

        const auto& si = contour[i];
        const std::size_t gr_idx = m_groups[k][i];

        if (gr_idx == std::size_t(-1)) {
          segments.push_back(si); continue;
        }

        const std::size_t im = i - 1;
        const std::size_t ip = i + 1;

        auto ss = si;
        const auto& sm = contour[im];
        const auto& sp = contour[ip];

        const auto& inm = m_saved[k][im];
        const auto& inp = m_saved[k][ip];

        correct_segment(sm, ss, sp, inm, inp);
        segments.push_back(ss);
      }
      segments.push_back(contour[n - 1]);
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

    FT get_angle_2(const FT angle) {

      FT angle_2 = angle;
      if (angle_2 > FT(90)) angle_2 = FT(180) - angle_2;
      else if (angle_2 < -FT(90)) angle_2 = FT(180) + angle_2;
      return angle_2;
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
  };

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_INTERNAL_SEGMENT_REGULARIZER_H
