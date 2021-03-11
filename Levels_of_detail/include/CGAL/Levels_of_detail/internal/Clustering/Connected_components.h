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

#ifndef CGAL_LEVELS_OF_DETAIL_INTERNAL_CLUSTERING_CONNECTED_COMPONENTS_H
#define CGAL_LEVELS_OF_DETAIL_INTERNAL_CLUSTERING_CONNECTED_COMPONENTS_H

#include <CGAL/license/Levels_of_detail.h>

// STL includes.
#include <vector>
#include <utility>

// Internal includes.
#include <CGAL/Levels_of_detail/internal/Spatial_search/Sphere_neighbor_query.h>
#include <CGAL/Levels_of_detail/internal/Shape_detection/Connected_component_region.h>
#include <CGAL/Levels_of_detail/internal/Shape_detection/Region_growing.h>

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

  template<
  typename GeomTraits,
  typename InputRange,
  typename PointMap>
  class Connected_components {

  public:
    using Traits = GeomTraits;
    using Input_range = InputRange;
    using Point_map = PointMap;

    using FT = typename Traits::FT;
    using Indices = std::vector<std::size_t>;

    using Neighbor_query =
    internal::Sphere_neighbor_query<Traits, Input_range, Point_map>;
    using Connected_component_region =
    internal::Connected_component_region;
    using Region_growing =
    internal::Region_growing<Input_range, Neighbor_query, Connected_component_region>;

    Connected_components(
      const Input_range& input_range,
      const Point_map& point_map,
      const FT scale,
      const std::size_t min_cluster_size) :
    m_input_range(input_range),
    m_point_map(point_map),
    m_scale(scale),
    m_min_cluster_size(min_cluster_size) {

      CGAL_precondition(m_input_range.size() > 0);
    }

    template<typename Cluster>
    void create_clusters(std::vector<Cluster>& clusters) const {

      Neighbor_query neighbor_query(
        m_input_range, m_scale, m_point_map);
      Connected_component_region region_type(m_min_cluster_size);
      Region_growing region_growing(
        m_input_range, neighbor_query, region_type);

      std::vector<Indices> regions;
      region_growing.detect(std::back_inserter(regions));

      clusters.clear();
      clusters.reserve(regions.size());

      Cluster cluster;
      for (const auto& region : regions) {
        cluster.clear();
        for (const std::size_t idx : region)
          cluster.push_back(*(m_input_range.begin() + idx));
        clusters.push_back(cluster);
      }
    }

  private:
    const Input_range& m_input_range;
    const Point_map& m_point_map;
    const FT m_scale;
    const std::size_t m_min_cluster_size;
  };

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_INTERNAL_CLUSTERING_CONNECTED_COMPONENTS_H
