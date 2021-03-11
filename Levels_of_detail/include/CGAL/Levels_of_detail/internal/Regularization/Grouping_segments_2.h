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

#ifndef CGAL_LEVELS_OF_DETAIL_INTERNAL_GROUPING_SEGMENTS_2
#define CGAL_LEVELS_OF_DETAIL_INTERNAL_GROUPING_SEGMENTS_2

#include <CGAL/license/Levels_of_detail.h>

#include <vector>
#include <map>
#include <utility>

#include <CGAL/Levels_of_detail/internal/utils.h>
#include <CGAL/Levels_of_detail/internal/Regularization/Segment_data_2.h>

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

  template<
    typename GeomTraits,
    typename Conditions>
  class Grouping_segments_2 {
  public:
    using Traits = GeomTraits;
    using FT = typename GeomTraits::FT;
    using Segment_data = typename internal::Segment_data_2<Traits>;
    using Targets_map = std::map <std::pair<std::size_t, std::size_t>, std::pair<FT, std::size_t>>;
    using Relations_map = std::map <std::pair<std::size_t, std::size_t>, std::pair<int, std::size_t>>;

    Grouping_segments_2() :
    m_tolerance(FT(1) / FT(1000000)),
    m_moe(FT(0)) { }

    void make_groups(const FT max_bound, const std::size_t n,
                       const std::map <std::size_t, Segment_data> & segments,
                       const std::vector<FT> & qp_result,
                       std::map<FT, std::vector<std::size_t>> & groups_by_value,
                       const Targets_map & targets, const Relations_map & relations = Relations_map()) {
      CGAL_precondition(n > 0);
      CGAL_precondition(max_bound > 0);
      CGAL_precondition(qp_result.size() > 0);

      m_cond.set_margin_of_error(max_bound);
      m_moe = m_cond.get_margin_of_error();
      CGAL_postcondition(m_moe > 0);

      groups_by_value.clear();
      m_groups.clear();
      m_segments_to_groups_hashmap.clear();
      m_values.clear();

      for (const auto & it : segments) {
        std::size_t seg_index = it.second.m_index;
        m_segments_to_groups_hashmap[seg_index] = -1;
      }

      build_initial_groups(n, targets, relations, qp_result);
      build_map_of_values(qp_result, segments);

      // Try to assign segments whose orientation has not been optimized thanks
      // to the regularization process, to an existing group.
      assign_segments_to_groups(segments);
      build_groups_by_value(groups_by_value);
    }

  private:
    const FT m_tolerance;
    FT m_moe;
    Conditions m_cond;
    std::map<std::size_t, int> m_segments_to_groups_hashmap;
    std::map <std::size_t, std::vector<std::size_t>> m_groups;
    std::map<int, FT> m_values;

    void build_initial_groups(const std::size_t n,
                              const Targets_map & targets, const Relations_map & relations,
                              const std::vector<FT> & qp_result) {
      std::size_t g = 0;
      auto rel_it = relations.begin();

      for (const auto & tar_it : targets) {
        const std::size_t i = tar_it.first.first;
        const std::size_t j = tar_it.first.second;
        const std::size_t p = tar_it.second.second;

        int r = 0;
        if (rel_it != relations.end()) {
          CGAL_precondition(rel_it->second.second == p);
          r = rel_it->second.first;
        }
        CGAL_postcondition(r == 0 || r == 1);

        if (CGAL::abs(qp_result[n + p]) >= m_tolerance) {
          if(rel_it != relations.end()) ++rel_it;
          continue;
        }

        const int g_i = m_segments_to_groups_hashmap[i];
        const int g_j = m_segments_to_groups_hashmap[j];
        const int groups_status = check_group_status(g_i, g_j);

        switch (groups_status) {
          case -1: break;

          case 1:
            r == 0 ? create_single_group(i, j, g) : create_separate_groups(i, j, g);
            break;

          case 2:
            r == 0 ? assign_segment_to_group(i, j) : create_new_group(i, g);
            break;

          case 3:
            r == 0 ? assign_segment_to_group(j, i) : create_new_group(j, g);
            break;

          case 4:
            if (r == 0) merge_two_groups(g_i, g_j);
            break;
        }

        if(rel_it != relations.end()) ++rel_it;
      }
    }

    void build_map_of_values(const std::vector<FT> & qp_result,
                             const std::map <std::size_t, Segment_data> & segments) {
      for (const auto & sm_i : m_segments_to_groups_hashmap) {
        int g_i = sm_i.second;

        if (g_i != -1 && (m_values.find(g_i) == m_values.end())) {
          const std::size_t seg_index = sm_i.first;
          const Segment_data & seg_data = segments.at(seg_index);
          const FT val = m_cond.reference(seg_data, qp_result[seg_index]);

          // Check if the angle that seems to be associated to this group of
          // segments is not too close to another value.
          int g_j = -1;
          for (const auto & it_m : m_values) {
            if (CGAL::abs(it_m.second - val) < m_moe)
              g_j = it_m.first;
          }

          if (g_j == -1)
            m_values[g_i] = val;
          else
            merge_two_groups(g_j, g_i);
        }
      }
    }

    void assign_segments_to_groups(const std::map <std::size_t, Segment_data> & segments) {
      for (const auto & sm_i : m_segments_to_groups_hashmap) {
        int g_i = sm_i.second;

        if (g_i == -1) {
          const std::size_t seg_index = sm_i.first;
          const Segment_data & seg_data = segments.at(seg_index);
          const FT val = m_cond.reference(seg_data, 0);
          int g_j = -1;

          for (const auto & it_m : m_values) {
            const FT val_j = it_m.second;
            const int g_index = it_m.first;

            g_j = m_cond.group_index(val, val_j, g_index);
            if (g_j != -1) break;
          }

          if (g_j == -1) {
            m_values.size() > 0 ? g_i = m_values.rbegin()->first + 1 : g_i = 0;
            m_values[g_i] = val;
          }
          else
            g_i = g_j;

          m_segments_to_groups_hashmap[seg_index] = g_i;
          m_groups[g_i].push_back(seg_index);
        }
      }
    }

    void build_groups_by_value(std::map <FT, std::vector<std::size_t>> & groups_by_value) {
      for (const auto & it_m : m_values) {
        const FT val = it_m.second;
        if (groups_by_value.find(val) == groups_by_value.end())
          groups_by_value[val] = std::vector<std::size_t>();
      }

      for (const auto & sm_i : m_segments_to_groups_hashmap) {
        const FT val = m_values.at(sm_i.second);
        if (groups_by_value.find(val) != groups_by_value.end())
          groups_by_value[val].push_back(sm_i.first);
      }
    }

    int check_group_status(const int g_i, const int g_j) const {
      if (g_i == -1 && g_j == -1) return 1;
      if (g_i == -1 && g_j != -1) return 2;
      if (g_i != -1 && g_j == -1) return 3;
      if (g_i != -1 && g_j != -1 && g_i != g_j) return 4;
      return -1;
    }

    void create_single_group (const std::size_t i, const std::size_t j, std::size_t & g) {
      m_segments_to_groups_hashmap[i] = g;
      m_segments_to_groups_hashmap[j] = g;
      m_groups[g].push_back(i);
      m_groups[g].push_back(j);
      ++g;
    }

    void create_separate_groups(const std::size_t i, const std::size_t j, std::size_t & g) {
      create_new_group(i, g);
      create_new_group(j, g);
    }

    void assign_segment_to_group(const std::size_t i, const std::size_t j) {
      const int g_j = m_segments_to_groups_hashmap[j];
      m_segments_to_groups_hashmap[i] = g_j;
      m_groups[g_j].push_back(i);
    }

    void create_new_group(const std::size_t i, std::size_t & g) {
      m_segments_to_groups_hashmap[i] = g;
      m_groups[g].push_back(i);
      ++g;
    }

    void merge_two_groups(const int g_i, const int g_j) {
      for (const auto gr : m_groups[g_j]) {
        m_segments_to_groups_hashmap[gr] = g_i;
        m_groups[g_i].push_back(gr);
      }
      m_groups[g_j].clear();
    }
  };

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_INTERNAL_GROUPING_SEGMENTS_2
