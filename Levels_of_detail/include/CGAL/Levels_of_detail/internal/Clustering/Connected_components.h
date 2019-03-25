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
#include <CGAL/Levels_of_detail/internal/Spacial_search/Sphere_neighbor_query.h>

#include <CGAL/Levels_of_detail/internal/Shape_detection/Estimate_normals_3.h>
#include <CGAL/Levels_of_detail/internal/Shape_detection/Connected_component_region.h>
#include <CGAL/Levels_of_detail/internal/Shape_detection/Region_growing.h>

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

  template<
  typename GeomTraits,
  typename ItemRange,
  typename PointMap>
  class Connected_components {

  public:
    using Traits = GeomTraits;
    using Item_range = ItemRange;
    using Point_map = PointMap;

    using FT = typename Traits::FT;

    using Neighbor_query = 
    internal::Sphere_neighbor_query<Traits, Item_range, Point_map>;
    using Normals_3 = 
    internal::Estimate_normals_3<Traits, Item_range, Point_map, Neighbor_query>;
    using Connected_component_region = internal::Connected_component_region;
    using Region_growing = 
    internal::Region_growing<Item_range, Neighbor_query, Connected_component_region>;

    Connected_components(
      const Item_range& items,
      const Point_map& point_map,
      const FT scale) : 
    m_items(items),
    m_point_map(point_map),
    m_scale(scale)
    { }

    template<typename OutputIterator>
    void create_clusters(OutputIterator clusters) const {
      
      Neighbor_query neighbor_query(
        m_items, m_scale, m_point_map);
      Normals_3 estimator(
        m_items, neighbor_query, m_point_map);
      Connected_component_region region(10);
      Region_growing region_growing(
        m_items, neighbor_query, region);
      region_growing.detect(clusters);
    }

  private:
    const Item_range& m_items;
    const Point_map& m_point_map;
    const FT m_scale;
  };

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_INTERNAL_CLUSTERING_CONNECTED_COMPONENTS_H
