// Copyright (c) 2019 INRIA Sophia-Antipolis (France).
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
//
// Author(s)     : Dmitry Anisimov, Simon Giraudot, Pierre Alliez, Florent Lafarge, and Andreas Fabri
//

#ifndef CGAL_LEVELS_OF_DETAIL_INTERNAL_SHAPE_DETECTION_REGION_GROWING_H
#define CGAL_LEVELS_OF_DETAIL_INTERNAL_SHAPE_DETECTION_REGION_GROWING_H

#include <CGAL/license/Levels_of_detail.h>

// STL includes.
#include <queue>
#include <vector>

// CGAL includes.
#include <CGAL/assertions.h>
#include <CGAL/property_map.h>

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

  template<
  typename InputRange,
  typename NeighborQuery,
  typename RegionType,
  typename SeedMap = CGAL::Identity_property_map<std::size_t> >
  class Region_growing {

  public:
    using Input_range = InputRange;
    using Neighbor_query = NeighborQuery;
    using Region_type = RegionType;
    using Seed_map = SeedMap;

    using Visited_items = std::vector<bool>;
    using Running_queue = std::queue<std::size_t>;
    using Indices       = std::vector<std::size_t>;

    Region_growing(
      const Input_range& input_range,
      Neighbor_query& neighbor_query,
      Region_type& region_type,
      const Seed_map seed_map = Seed_map()) :
    m_input_range(input_range),
    m_neighbor_query(neighbor_query),
    m_region_type(region_type),
    m_seed_map(seed_map),
    m_use_overlap(false) {

      CGAL_precondition(m_input_range.size() > 0);
      clear();
    }

    template<typename OutputIterator>
    OutputIterator detect(OutputIterator regions) {

      clear();
      Indices region;

      // Grow regions.
      for (std::size_t i = 0; i < m_input_range.size(); ++i) {
        const std::size_t seed_index = get(m_seed_map, i);

        // Skip items that user does not want to use.
        if (seed_index == std::size_t(-1))
          continue;

        CGAL_precondition(
          seed_index >= 0 && seed_index < m_input_range.size());

        // Try to grow a new region from the index of the seed item.
        if (!m_visited[seed_index]) {
          propagate(seed_index, region);

          // Check global conditions.
          if (!m_region_type.is_valid_region(region))
            revert(region);
          else
            *(regions++) = region;
        }
      }
      return regions;
    }

    template<typename OutputIterator>
    OutputIterator unassigned_items(OutputIterator output) const {

      // Return indices of all unassigned items.
      for (std::size_t i = 0; i < m_input_range.size(); ++i) {
        const std::size_t seed_index = get(m_seed_map, i);

        // Skip items that user does not want to use.
        if (seed_index == std::size_t(-1))
          continue;

        CGAL_precondition(
          seed_index >= 0 && seed_index < m_input_range.size());

        if (!m_visited[seed_index])
          *(output++) = seed_index;
      }
      return output;
    }

    void clear() {

      m_visited.clear();
      m_visited.resize(m_input_range.size(), false);

      if (m_use_overlap) {
        m_extra.clear();
        m_extra.resize(m_input_range.size(), false);
      }
    }

    void release_memory() {

      m_visited.clear();
      m_visited.shrink_to_fit();
    }

    void use_overlap(const bool state) {
      m_use_overlap = state;
    }

  private:
    bool m_use_overlap;

    void propagate(const std::size_t seed_index, Indices& region) {
      region.clear();

      // Use two queues, while running on this queue, push to the other queue;
      // When the queue is done, update the shape of the current region and swap to the other queue;
      // depth_index is the index of the queue we are using.
      Running_queue running_queue[2];
      bool depth_index = 0;

      // Once the index of an item is pushed to the queue, it is pushed to the region too.
      m_visited[seed_index] = true;
      running_queue[depth_index].push(seed_index);
      region.push_back(seed_index);

      // Update internal properties of the region.
      m_region_type.update(region);

      if (m_use_overlap) {
        m_extra.clear();
        m_extra.resize(m_input_range.size(), false);
        m_extra[seed_index] = true;
      }

      Indices neighbors;
      while (
        !running_queue[depth_index].empty() ||
        !running_queue[!depth_index].empty()) {

        // Call the next item index of the queue and remove it from the queue.
        const std::size_t item_index = running_queue[depth_index].front();
        running_queue[depth_index].pop();

        // Get neighbors of the current item.
        neighbors.clear();
        m_neighbor_query(item_index, neighbors);

        // Visit all found neighbors.
        for (const std::size_t neighbor_index : neighbors) {

          // Skip items that user does not want to use.
          if (neighbor_index == std::size_t(-1))
            continue;

          CGAL_precondition(
            neighbor_index >= 0 && neighbor_index < m_input_range.size());

          if (!m_use_overlap) {

            if (!m_visited[neighbor_index] &&
              m_region_type.is_part_of_region(item_index, neighbor_index, region)) {

              // Add this neighbor to the other queue so that we can visit it later.
              m_visited[neighbor_index] = true;
              running_queue[!depth_index].push(neighbor_index);
              region.push_back(neighbor_index);
            }

          } else {

            if (m_region_type.is_already_visited(item_index, neighbor_index, m_visited[neighbor_index])) {

              m_visited[neighbor_index] = true;
              running_queue[!depth_index].push(neighbor_index);
              region.push_back(neighbor_index);

            } else {

              if (!m_extra[neighbor_index] &&
              m_region_type.is_part_of_region(item_index, neighbor_index, region)) {

                running_queue[!depth_index].push(neighbor_index);
                m_extra[neighbor_index] = true;
              }
            }
          }
        }

        // Update internal properties of the region.
        if (running_queue[depth_index].empty()) {

          m_region_type.update(region);
          depth_index = !depth_index;
        }
      }
    }

    void revert(const Indices& region) {
      for (const std::size_t item_index : region)
        m_visited[item_index] = false;
    }

    const Input_range& m_input_range;
    Neighbor_query& m_neighbor_query;
    Region_type& m_region_type;
    const Seed_map m_seed_map;

    Visited_items m_visited;
    Visited_items m_extra;
  };

} // namespace internal
} // namespace Levels_of_detail
} // namespace CGAL

#endif // CGAL_LEVELS_OF_DETAIL_INTERNAL_SHAPE_DETECTION_REGION_GROWING_H
