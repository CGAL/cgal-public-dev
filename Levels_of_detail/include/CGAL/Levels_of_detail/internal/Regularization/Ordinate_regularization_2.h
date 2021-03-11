// Copyright (c) 2019 GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
//
// Author(s)     : Jean-Philippe Bauchet, Florent Lafarge, Gennadii Sytov, Dmitry Anisimov
//

#ifndef CGAL_LEVELS_OF_DETAIL_INTERNAL_ORDINATE_REGULARIZATION_2
#define CGAL_LEVELS_OF_DETAIL_INTERNAL_ORDINATE_REGULARIZATION_2

#include <CGAL/license/Levels_of_detail.h>

#include <map>
#include <utility>
#include <vector>
#include <iostream>

#include <CGAL/Levels_of_detail/internal/utils.h>
#include <CGAL/Levels_of_detail/internal/Regularization/Segment_data_2.h>
#include <CGAL/Levels_of_detail/internal/Regularization/Grouping_segments_2.h>
#include <CGAL/Levels_of_detail/internal/Regularization/Conditions_ordinates_2.h>

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

  template<
    typename GeomTraits,
    typename InputRange,
    typename SegmentMap>
  class Ordinate_regularization_2 {
  public:

    using Traits = GeomTraits;
    using Input_range = InputRange;
    using Segment_map = SegmentMap;

    typedef typename GeomTraits::FT FT;

    using Segment = typename GeomTraits::Segment_2;
    using Point = typename GeomTraits::Point_2;
    using Segment_data = typename internal::Segment_data_2<Traits>;
    using Conditions = typename internal::Conditions_ordinates_2<Traits>;
    using Grouping = internal::Grouping_segments_2<Traits, Conditions>;
    using Vector  = typename GeomTraits::Vector_2;
    using Targets_map = std::map <std::pair<std::size_t, std::size_t>, std::pair<FT, std::size_t>>;

    Ordinate_regularization_2 (
      InputRange& input_range,
      const FT d_max = FT(0.1),
      const SegmentMap segment_map = SegmentMap()) :
    m_input_range(input_range),
    m_d_max(CGAL::abs(d_max)),
    m_segment_map(segment_map),
    m_modified_segments_counter(0) { }

    FT target_value(const std::size_t i, const std::size_t j) {
      if(m_segments.size() == 0) return FT(0);
      CGAL_precondition(m_segments.size() > 0);
      CGAL_precondition(m_segments.find(i) != m_segments.end());
      CGAL_precondition(m_segments.find(j) != m_segments.end());

      const Segment_data & s_i = m_segments.at(i);
      const Segment_data & s_j = m_segments.at(j);

      const FT tar_val = s_i.m_reference_coordinates.y() - s_j.m_reference_coordinates.y();
      if (CGAL::abs(tar_val) < bound(i) + bound(j)) {
        m_targets[std::make_pair(i, j)] = tar_val;
      }

      return tar_val;
    }

    FT bound(const std::size_t i) const {
      CGAL_precondition(i >= 0 && i < m_input_range.size());
      return m_d_max;
    }

    void update(const std::vector<FT> & result) {
      CGAL_precondition(result.size() > 0);
      const std::size_t n = m_input_range.size();
      std::map <FT, std::vector<std::size_t>> collinear_groups_by_ordinates;
      std::map <std::size_t, Segment_data> segments;
      Targets_map targets;

      for (const auto & group : m_parallel_groups) {
        if (group.size() < 2) continue;

        collinear_groups_by_ordinates.clear();
        segments.clear();
        targets.clear();

        build_grouping_data(group, segments, targets);

        if (segments.size() > 0) {
          m_grouping.make_groups(m_d_max, n, segments, result, collinear_groups_by_ordinates, targets);
          translate_collinear_segments(collinear_groups_by_ordinates);
        }
      }
    }

    template<typename Range, typename IndexMap = CGAL::Identity_property_map<std::size_t>>
  	void add_group(const Range& group, const IndexMap index_map = IndexMap()) {
      std::vector<std::size_t> gr;
      for (const auto & item : group) {
        const std::size_t seg_index = get(index_map, item);
        gr.push_back(seg_index);
      }

      if (gr.size() > 1) {
        m_parallel_groups.push_back(gr);
        build_segment_data_map(gr);
      }
    }

    std::size_t number_of_modified_segments() const {
      return m_modified_segments_counter;
    }

    void make_bounds() {
      m_bounds.clear();
      m_bounds.resize(m_input_range.size(), m_d_max);
      std::vector<FT> norms(m_input_range.size(), FT(0));

      FT max_length = -FT(1);
      for (const auto& segment : m_segments) {
        const FT length = segment.second.m_length;
        max_length = CGAL::max(max_length, length);
        norms[segment.second.m_index] = length;
      }

      for (std::size_t i = 0; i < norms.size(); ++i) {
        norms[i] /= max_length;
        norms[i] *= m_d_max;
        m_bounds[i] -= norms[i];
      }
    }

  private:
    Input_range& m_input_range;
    const FT m_d_max;
    const Segment_map  m_segment_map;
    std::map <std::size_t, Segment_data> m_segments;
    std::map <std::pair<std::size_t, std::size_t>, FT> m_targets;
    Grouping m_grouping;
    std::vector <std::vector <std::size_t>> m_parallel_groups;
    std::size_t m_modified_segments_counter;
    std::vector<FT> m_bounds;

    void build_segment_data_map(const std::vector<std::size_t> & paral_gr) {
      if (paral_gr.size() < 2) return;

      Point frame_origin;
      for(std::size_t i = 0; i < paral_gr.size(); ++i) {
        const std::size_t seg_index = paral_gr[i];

        CGAL_precondition(m_segments.find(seg_index) == m_segments.end());
        if(m_segments.find(seg_index) != m_segments.end())
          continue;

        const Segment& seg = get(m_segment_map, *(m_input_range.begin() + seg_index));
        Segment_data seg_data(seg, seg_index);

        if (i == 0)
          frame_origin = seg_data.m_barycentre;

        seg_data.m_reference_coordinates = internal::transform_coordinates(
                seg_data.m_barycentre, frame_origin, seg_data.m_orientation);
        m_segments.emplace(seg_index, seg_data);
      }
    }

    void build_grouping_data(const std::vector <std::size_t> & group,
                             std::map <std::size_t, Segment_data> & segments,
                             Targets_map & targets) {
      for (const std::size_t it : group) {
        const std::size_t seg_index = it;

        CGAL_precondition(m_segments.find(seg_index) != m_segments.end());
        const Segment_data& seg_data = m_segments.at(seg_index);

        segments.emplace(seg_index, seg_data);
        std::size_t tar_index = 0;

        for(const auto & ti : m_targets) {
          const std::size_t seg_index_tar_i = ti.first.first;
          const std::size_t seg_index_tar_j = ti.first.second;
          const FT tar_val = ti.second;

          if (seg_index_tar_i == seg_index) {
            targets[std::make_pair(seg_index_tar_i, seg_index_tar_j)] = std::make_pair(tar_val, tar_index);
          }

          ++tar_index;
        }
      }
    }

    void translate_collinear_segments(const std::map <FT, std::vector<std::size_t>> & collinear_groups_by_ordinates) {
      for (const auto & mi : collinear_groups_by_ordinates) {
        const FT dt = mi.first;
        const std::vector<std::size_t> & group = mi.second;
        int l_index = find_longest_segment(group);
        CGAL_postcondition(l_index >= 0);
        if(l_index < 0) {
          std::cerr << "Cannot translate collinear segments! Cannot find the longest segment!" << std::endl;
          return;
        }

        CGAL_precondition(m_segments.find(l_index) != m_segments.end());
        const Segment_data & l_data = m_segments.at(l_index);

        FT new_difference = dt - l_data.m_reference_coordinates.y();
        set_difference(l_index, new_difference);

        const FT l_a = l_data.m_a;
        const FT l_b = l_data.m_b;
        const FT l_c = l_data.m_c;
        const Vector & l_direction = l_data.m_direction;

        // Translate the other segments, so that they rest upon the line ax + by + c = 0.
        for (const std::size_t it : group) {
          if (it != static_cast<std::size_t>(l_index)) {
            CGAL_precondition(m_segments.find(it) != m_segments.end());
            const Segment_data & seg_data = m_segments.at(it);

            new_difference = dt - seg_data.m_reference_coordinates.y();
            set_difference(it, new_difference, l_a, l_b, l_c, l_direction);
          }
        }
      }
    }

    int find_longest_segment(const std::vector<std::size_t> & group) const {
      FT l_max = -FT(1000000000000);
      int l_index = -1;

      for (const std::size_t it : group) {
        const FT seg_length = m_segments.at(it).m_length;

        if (l_max < seg_length) {
          l_max = seg_length;
          l_index = it;
        }
      }
      return l_index;
    }

    void set_difference(const int i, const FT new_difference) {
      const FT difference = new_difference;
      Segment_data & seg_data = m_segments.at(i);

      const Vector & direction = seg_data.m_direction;
      const Vector final_normal = Vector(-direction.y(), direction.x());

      const Point &source = seg_data.m_segment.source();
      const Point &target = seg_data.m_segment.target();

      Point new_source = Point(source.x() + difference * final_normal.x(), source.y() + difference * final_normal.y());
      Point new_target = Point(target.x() + difference * final_normal.x(), target.y() + difference * final_normal.y());

      const FT bx = (new_source.x() + new_target.x()) / FT(2);
      const FT by = (new_source.y() + new_target.y()) / FT(2);

      m_input_range[i] = Segment(new_source, new_target);
      seg_data.m_c = -seg_data.m_a * bx - seg_data.m_b * by;

      ++m_modified_segments_counter;
    }

    void set_difference(const int i, const FT new_difference, const FT a, const FT b, const FT c, const Vector &direction) {
      FT difference = new_difference;
      Segment_data & seg_data = m_segments.at(i);

      seg_data.m_direction = direction;
      if (seg_data.m_direction.y() < FT(0) || (seg_data.m_direction.y() == FT(0) && seg_data.m_direction.x() < FT(0)))
        seg_data.m_direction = -seg_data.m_direction;

      Vector final_normal = Vector(-seg_data.m_direction.y(), seg_data.m_direction.x());
      FT x1, x2, y1, y2;

      const Point &source = seg_data.m_segment.source();
      const Point &target = seg_data.m_segment.target();

      if (CGAL::abs(seg_data.m_direction.x()) > CGAL::abs(seg_data.m_direction.y())) {
        x1 = source.x() + difference * final_normal.x();
        x2 = target.x() + difference * final_normal.x();

        y1 = (-c - a * x1) / b;
        y2 = (-c - a * x2) / b;
      }
      else {
        y1 = source.y() + difference * final_normal.y();
        y2 = target.y() + difference * final_normal.y();

        x1 = (-c - b * y1) / a;
        x2 = (-c - b * y2) / a;
      }

      const Point new_source = Point(x1, y1);
      const Point new_target = Point(x2, y2);
      m_input_range[i] = Segment(new_source, new_target);

      ++m_modified_segments_counter;
    }
  };

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_INTERNAL_ORDINATE_REGULARIZATION_2
