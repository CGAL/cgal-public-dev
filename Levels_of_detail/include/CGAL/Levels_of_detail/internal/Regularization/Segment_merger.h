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

#ifndef CGAL_LEVELS_OF_DETAIL_INTERNAL_SEGMENT_MERGER_H
#define CGAL_LEVELS_OF_DETAIL_INTERNAL_SEGMENT_MERGER_H

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
  class Segment_merger {

  public:
    using Traits = GeomTraits;

    using FT = typename Traits::FT;
    using Point_2 = typename Traits::Point_2;
    using Point_3 = typename Traits::Point_3;
    using Segment_2 = typename Traits::Segment_2;
    using Vector_2 = typename Traits::Vector_2;
    using Line_2 = typename Traits::Line_2;

    using Indices = std::vector<std::size_t>;
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

    Segment_merger(
      const FT ordinate_bound,
      const FT angle_threshold = FT(5)) :
    m_ordinate_bound(ordinate_bound),
    m_angle_threshold(angle_threshold),
    m_pi(static_cast<FT>(CGAL_PI)),
    m_k(10),
    m_num_samples_per_segment(m_k * 2)
    { }

    void merge_segments_with_outer_boundary(
      const std::vector<Segment_2>& outer,
      std::vector<Segment_2>& segments) {

      auto merged = outer;
      const FT default_bound = FT(1) / FT(1000);
      merge_with_outer_boundary(default_bound, outer, merged);
      merge_with_outer_boundary(m_ordinate_bound, outer, segments);
      for (const auto& segment : segments)
        merged.push_back(segment);
      segments = merged;
    }

    void merge_segments(
      std::vector<Segment_2>& segments,
      const bool use_weighted = true) {

      std::vector<Segment_2> merged;

      std::size_t num_groups = 0;
      std::vector<Segment_2> group;
      std::vector<bool> states(segments.size(), false);

      for (std::size_t i = 0; i < segments.size(); ++i) {
        const auto& segment_i = segments[i];
        if (states[i]) continue;

        create_collinear_group(
          m_ordinate_bound, segments, segment_i, i, states, group);

        Segment_2 ref_segment;
        if (use_weighted) ref_segment = find_weighted_segment(group);
        else {
          const std::size_t seg_idx = find_longest_segment(group);
          ref_segment = group[seg_idx];
        }
        const bool success = create_merged_segment(group, ref_segment);

        if (success) {
          merged.push_back(ref_segment); ++num_groups;
        }
      }
      std::cout << "Num collinear groups (wrt inner): " << num_groups << std::endl;

      std::vector<Segment_2> clean;
      remove_zero_length_segments(merged, clean);
      segments = clean;
    }

    void snap_segments(
      const std::vector<Segment_2>& segments_outer,
      std::vector<Segment_2>& segments_inner) {

      merge_with_outer_boundary(
        m_ordinate_bound, segments_outer, segments_inner);
      connect_to_corners(
        segments_outer, segments_inner);

      std::vector<Segment_2> clean;
      remove_zero_length_segments(segments_inner, clean);
      segments_inner = clean;
    }

  private:
    const FT m_ordinate_bound;
    const FT m_angle_threshold;
    const FT m_pi;

    const std::size_t m_k;
    const std::size_t m_num_samples_per_segment;

    void merge_with_outer_boundary(
      const FT ordinate_bound,
      const std::vector<Segment_2>& segments_outer,
      std::vector<Segment_2>& segments_inner) {

      std::vector<Segment_2> merged;

      std::vector<bool> states(segments_inner.size(), false);
      std::vector<Segment_2> group;
      std::size_t num_groups = 0;

      for (std::size_t i = 0; i < segments_outer.size(); ++i) {
        const auto& segment_i = segments_outer[i];

        create_collinear_group(
          ordinate_bound,
          segments_inner, segment_i, std::size_t(-1), states, group);

        if (group.size() > 0) {
          Segment_2 ref_segment = segment_i;
          const bool success = create_merged_segment(group, ref_segment);

          if (success) {
            merged.push_back(ref_segment);
            ++num_groups;
          }
        }
      }

      for (std::size_t i = 0; i < segments_inner.size(); ++i) {
        if (states[i]) continue;
        merged.push_back(segments_inner[i]);
      }

      segments_inner = merged;
      std::cout << "Num collinear groups (wrt outer): " << num_groups << std::endl;
    }

    void create_collinear_group(
      const FT ordinate_bound,
      const std::vector<Segment_2>& segments,
      const Segment_2& segment_i,
      const std::size_t i,
      std::vector<bool>& states,
      std::vector<Segment_2>& group) {

      group.clear();
      if (i != std::size_t(-1)) {

        group.push_back(segment_i);
        states[i] = true;
      }

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

          if (distance <= ordinate_bound) {
            group.push_back(segment_j); states[j] = true;
          }
        }
      }
    }

    bool create_merged_segment(
      const std::vector<Segment_2>& group,
      Segment_2& ref_segment) {

      const Line_2 line =
        Line_2(ref_segment.source(), ref_segment.target());

      if (line.a() == FT(0) && line.b() == FT(0) && line.c() == FT(0))
        return false;

      std::vector<Point_2> points;
      for (const auto& segment : group) {

        const Point_2 p = line.projection(segment.source());
        const Point_2 q = line.projection(segment.target());

        points.push_back(p);
        points.push_back(q);
      }
      update_segment(points, ref_segment);
      return true;
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

      if (source == target) {
        const std::size_t seg_idx = find_longest_segment(segments);
        return segments[seg_idx];
      }
      return Segment_2(source, target);
    }

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

    void remove_zero_length_segments(
      const std::vector<Segment_2>& contour,
      std::vector<Segment_2>& segments) {

      segments.clear();
      for (const auto& segment : contour)
        if (segment.squared_length() > internal::tolerance<FT>())
          segments.push_back(segment);
    }

    void connect_to_corners(
      const std::vector<Segment_2>& segments_outer,
      std::vector<Segment_2>& segments_inner) {

      std::vector<Point_pair> pair_range;
      create_pair_range(segments_outer, pair_range);

      Point_map point_map;
      K_neighbor_query neighbor_query(pair_range, FT(m_k), point_map);

      Indices neighbors;
      for (auto& segment : segments_inner) {
        auto s = segment.source();
        auto t = segment.target();

        neighbor_query(s, neighbors);
        update_point(
          segments_outer, pair_range, neighbors, s);

        neighbor_query(t, neighbors);
        update_point(
          segments_outer, pair_range, neighbors, t);

        segment = Segment_2(s, t);
      }
    }

    void create_pair_range(
      const std::vector<Segment_2>& segments_outer,
      std::vector<Point_pair>& pair_range) {

      create_range(segments_outer, pair_range);
      find_corners(segments_outer, pair_range);
    }

    void create_range(
      const std::vector<Segment_2>& segments_outer,
      std::vector<Point_pair>& pair_range) {

      pair_range.clear();
      std::vector<Point_2> samples; Range_data data;
      for (std::size_t i = 0; i < segments_outer.size(); ++i) {
        const auto& segment = segments_outer[i];

        const auto& s = segment.source();
        const auto& t = segment.target();

        samples.clear();
        Point_generator generator(s, t, m_num_samples_per_segment);
        std::copy_n(generator, m_num_samples_per_segment - 1,
        std::back_inserter(samples));

        data.is_corner = true;
        data.seg_i = i;
        pair_range.push_back(std::make_pair(samples[0], data));

        for (std::size_t j = 1; j < samples.size(); ++j) {
          data.is_corner = false;
          data.seg_i = i;
          pair_range.push_back(std::make_pair(samples[j], data));
        }
      }
    }

    void find_corners(
      const std::vector<Segment_2>& segments_outer,
      std::vector<Point_pair>& pair_range) {

      std::size_t num_corners = 0;
      for (auto& pair : pair_range) {

        const auto& p = pair.first;
        auto& data = pair.second;

        if (data.is_corner) {
          const bool success = find_corner(p, segments_outer, data);
          if (success) ++num_corners;
        }
      }
      std::cout << "Num corners: " << num_corners << std::endl;
    }

    bool find_corner(
      const Point_2& p,
      const std::vector<Segment_2>& segments_outer,
      Range_data& data) {

      const std::size_t n = segments_outer.size();
      for (std::size_t i = 0; i < n; ++i) {
        const auto& segment = segments_outer[i];
        const auto& s = segment.source();

        if (p == s) {
          const std::size_t im = (i + n - 1) % n;
          data.seg_i = im;
          data.seg_j = i;
          return true;
        }
      }
      return false;
    }

    void update_point(
      const std::vector<Segment_2>& segments_outer,
      const std::vector<Point_pair>& pair_range,
      const Indices& neighbors,
      Point_2& query) {

      for (const std::size_t idx : neighbors) {

        const auto& p    = pair_range[idx].first;
        const auto& data = pair_range[idx].second;

        if (data.seg_i != std::size_t(-1) && data.seg_j != std::size_t(-1)) {

          const auto& segment_i = segments_outer[data.seg_i];
          const auto& segment_j = segments_outer[data.seg_j];

          const FT angle   = angle_degree_2(segment_i, segment_j);
          const FT angle_2 = get_angle_2(angle);

          if (
            CGAL::abs(angle_2) >= m_angle_threshold &&
            CGAL::abs(angle_2) <= FT(90) - m_angle_threshold) {

            if (internal::distance(p, query) <= m_ordinate_bound * FT(2)) {
              query = p; return;
            }
          }
        }
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

#endif // CGAL_LEVELS_OF_DETAIL_INTERNAL_SEGMENT_MERGER_H
