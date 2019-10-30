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

#ifndef CGAL_LEVELS_OF_DETAIL_INTERNAL_PARTITION_23_ADAPTER_H
#define CGAL_LEVELS_OF_DETAIL_INTERNAL_PARTITION_23_ADAPTER_H

// STL includes.
#include <set>
#include <map>
#include <vector>

// CGAL includes.
#include <CGAL/property_map.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/Random.h>

// Internal includes.
#include <CGAL/Levels_of_detail/internal/struct.h>
#include <CGAL/Levels_of_detail/internal/utils.h>

// Testing.
#include "../../../../../test/Levels_of_detail/include/Saver.h"

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

template<typename GeomTraits>
class Partition_23_adapter {

public:
  using Traits = GeomTraits;

  using FT = typename Traits::FT;
  using Point_2 = typename Traits::Point_2;
  using Point_3 = typename Traits::Point_3;
  using Segment_2 = typename Traits::Segment_2;
  using Plane_3 = typename Traits::Plane_3;
  using Triangle_3 = typename Traits::Triangle_3;

  using Indices = std::vector<std::size_t>;
  using Partition_2 = internal::Partition_2<Traits>;
  using Face = typename Partition_2::Face;
  using Random = CGAL::Random;
  using Point_generator = CGAL::Random_points_in_triangle_3<Point_3>;
  using Approximate_face = internal::Partition_edge_3<Traits>;

  Partition_23_adapter(
    const std::map<std::size_t, Plane_3>& plane_map,
    Partition_2& partition_2) :
  m_plane_map(plane_map),
  m_partition_2(partition_2),
  m_num_samples(100),
  m_random(0) 
  { }

  void get_points_for_visibility_3(
    std::vector<Point_3>& points,
    std::vector<Indices>& updated_regions,
    std::vector<Plane_3>& planes) {

    points.clear();
    updated_regions.clear();
    planes.clear();

    const std::size_t num_labels = update_labels();
    updated_regions.resize(num_labels);
    planes.resize(num_labels);

    std::size_t count = 0;
    for (const auto& pface : m_partition_2.faces) {
      if (pface.visibility == Visibility_label::OUTSIDE)
        continue;
      add_samples(pface, points, updated_regions, planes, count);
    }

    Saver<Traits> saver;
    saver.export_points(
      points, 
      updated_regions,
      "/Users/monet/Documents/lod/logs/buildings/tmp/visibility_points_3");
  }

  bool get_approximate_roofs(
    std::vector<Approximate_face>& building_roofs) {

    building_roofs.clear();
    return false;
  }

  bool get_approximate_walls(
    std::vector<Approximate_face>& building_walls) {

    building_walls.clear();
    return false;
  }

private:
  const std::map<std::size_t, Plane_3>& m_plane_map;
  Partition_2& m_partition_2;
  
  const std::size_t m_num_samples;
  Random m_random;

  std::size_t update_labels() {

    std::set<std::size_t> tmp;
    for (const auto& pface : m_partition_2.faces) {
      if (pface.visibility == Visibility_label::OUTSIDE)
        continue;
      tmp.insert(pface.label);
    }

    Indices labels;
    labels.reserve(tmp.size());

    for (const std::size_t label : tmp)
      labels.push_back(label);

    for (auto& pface : m_partition_2.faces) {
      if (pface.visibility == Visibility_label::OUTSIDE)
        continue;
      
      for (std::size_t i = 0; i < labels.size(); ++i) {
        if (pface.label == labels[i]) {
          pface.plane = m_plane_map.at(pface.label);
          pface.label = i; break;
        }
      }
    }
    return labels.size();
  }

  void add_samples(
    const Face& pface,
    std::vector<Point_3>& points,
    std::vector<Indices>& updated_regions,
    std::vector<Plane_3>& planes,
    std::size_t& count) {

    const auto& polygon = pface.outer_polygon;
    const auto& plane = pface.plane;
    const std::size_t label = pface.label;

    planes[label] = plane;

    const auto& p0 = polygon[0];
    const Point_3 q0 = internal::position_on_plane_3(p0, plane);

    std::vector<Point_3> samples;
    for (std::size_t i = 1; i < polygon.size() - 1; ++i) {
      const std::size_t ip = i + 1;

      const auto& p1 = polygon[i];
      const auto& p2 = polygon[ip];

      const Point_3 q1 = internal::position_on_plane_3(p1, plane);
      const Point_3 q2 = internal::position_on_plane_3(p2, plane);

      const Triangle_3 triangle = Triangle_3(q0, q1, q2);

      samples.clear();
      Point_generator generator(triangle, m_random);
      std::copy_n(
        generator, m_num_samples, std::back_inserter(samples));

      for (const auto& sample : samples) {
        points.push_back(sample);
        updated_regions[label].push_back(count);
        ++count;
      }
    }
  }
};

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_INTERNAL_PARTITION_23_ADAPTER_H
