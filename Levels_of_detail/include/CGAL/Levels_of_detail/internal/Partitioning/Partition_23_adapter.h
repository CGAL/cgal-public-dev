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
#include <CGAL/Random.h>
#include <CGAL/property_map.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/point_generators_3.h>

// Buildings.
#include <CGAL/Levels_of_detail/internal/Buildings/Building_roofs_estimator.h>
#include <CGAL/Levels_of_detail/internal/Buildings/Building_walls_estimator.h>

// Regularization.
#include <CGAL/Levels_of_detail/internal/Regularization/Segment_merger.h>

// Internal includes.
#include <CGAL/Levels_of_detail/internal/utils.h>
#include <CGAL/Levels_of_detail/internal/struct.h>

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

  using Random = CGAL::Random;
  using Points_3 = std::vector<Point_3>;
  using Indices = std::vector<std::size_t>;
  using Partition_2 = internal::Partition_2<Traits>;
  using Face = typename Partition_2::Face;
  using Approximate_face = internal::Partition_edge_3<Traits>;
  using Identity_map_3 = CGAL::Identity_property_map<Point_3>;
  using Boundary = internal::Boundary<Traits>;

  using Point_generator_seg = CGAL::Points_on_segment_2<Point_2>;
  using Point_generator_tri = CGAL::Random_points_in_triangle_3<Point_3>;

  using Building_roofs_estimator =
    internal::Building_roofs_estimator<Traits, Points_3, Identity_map_3>;
  using Building_walls_estimator =
    internal::Building_walls_estimator<Traits>;

  using Segment_merger = internal::Segment_merger<Traits>;

  Partition_23_adapter(
    const std::map<std::size_t, Plane_3>& plane_map,
    const FT bottom_z, const FT top_z,
    const FT ordinate_bound_2,
    const FT max_height_difference,
    Partition_2& partition_2) :
  m_plane_map(plane_map),
  m_bottom_z(bottom_z), m_top_z(top_z),
  m_ordinate_bound_2(ordinate_bound_2),
  m_max_height_difference(max_height_difference),
  m_partition_2(partition_2),
  m_num_samples_per_segment(10),
  m_num_samples_per_triangle(100),
  m_random(0) {
    m_num_labels = update_labels();
  }

  void get_points_for_visibility_3(
    std::vector<Point_3>& points,
    std::vector<Indices>& updated_regions,
    std::vector<Plane_3>& planes) {

    points.clear();
    updated_regions.clear();
    planes.clear();

    updated_regions.resize(m_num_labels);
    planes.resize(m_num_labels);

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
      "/Users/monet/Documents/gf/lod/logs/buildings/tmp/visibility_points_3");
  }

  bool get_approximate_roofs(
    std::vector<Approximate_face>& building_roofs) {

    building_roofs.clear();

    Points_3 points;
    std::vector<Indices> regions;
    Identity_map_3 identity_map_3;
    create_roof_data(points, regions);

    Building_roofs_estimator estimator(points, identity_map_3, regions);
    estimator.estimate(building_roofs);
    return true;
  }

  bool get_approximate_inner_walls(
    std::vector<Approximate_face>& building_inner_walls) {

    std::vector<Segment_2> segments;
    create_inner_segments(segments);
    merge_segments(segments);
    create_walls(segments, building_inner_walls);

    Saver<Traits> saver;
    saver.save_polylines(
      segments,
      "/Users/monet/Documents/gf/lod/logs/buildings/tmp/clean_internal_segments");

    return true;
  }

private:
  const std::map<std::size_t, Plane_3>& m_plane_map;
  const FT m_bottom_z, m_top_z;
  const FT m_ordinate_bound_2;
  const FT m_max_height_difference;
  Partition_2& m_partition_2;

  std::size_t m_num_labels;
  const std::size_t m_num_samples_per_segment;
  const std::size_t m_num_samples_per_triangle;
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
      Point_generator_tri generator(triangle, m_random);
      std::copy_n(
        generator, m_num_samples_per_triangle, std::back_inserter(samples));

      for (const auto& sample : samples) {
        points.push_back(sample);
        updated_regions[label].push_back(count);
        ++count;
      }
    }
  }

  void create_roof_data(
    Points_3& points,
    std::vector<Indices>& regions) {

    points.clear();

    regions.clear();
    regions.resize(m_num_labels);

    std::size_t count = 0;
    for (const auto& pface : m_partition_2.faces) {
      if (pface.visibility == Visibility_label::OUTSIDE)
        continue;

      const auto& plane = pface.plane;
      const std::size_t label = pface.label;
      const auto& polygon = pface.outer_polygon;

      for (const auto& p : polygon) {
        const Point_3 q = internal::position_on_plane_3(p, plane);
        points.push_back(q); regions[label].push_back(count); ++count;
      }
    }
  }

  void create_inner_segments(
    std::vector<Segment_2>& segments) {

    segments.clear();
    for (const auto& pface : m_partition_2.faces) {
      if (pface.visibility == Visibility_label::OUTSIDE)
        continue;

      const auto& neighbors = pface.neighbors;
      const auto& edges     = pface.edges;

      for (std::size_t i = 0; i < neighbors.size(); ++i) {
        const int idx = neighbors[i];
        if (idx < 0) continue;

        const auto& nface = m_partition_2.faces[idx];
        if (nface.visibility == Visibility_label::OUTSIDE)
          continue;
        if (nface.label == pface.label)
          continue;

        const auto& segment = edges[i];
        if (is_inner_edge(segment, pface.plane, nface.plane))
          segments.push_back(segment);
      }
    }
  }

  bool is_inner_edge(
    const Segment_2& segment,
    const Plane_3& plane1, const Plane_3& plane2) {

    const auto& s = segment.source();
    const auto& t = segment.target();

    std::vector<Point_2> samples;
    Point_generator_seg generator(s, t, m_num_samples_per_segment);
    std::copy_n(generator, m_num_samples_per_segment - 1,
    std::back_inserter(samples));

    FT max_diff = -FT(1);
    for (const auto& p : samples) {
      const Point_3 q1 = internal::position_on_plane_3(p, plane1);
      const Point_3 q2 = internal::position_on_plane_3(p, plane2);

      const FT diff = CGAL::abs(q1.z() - q2.z());
      max_diff = CGAL::max(diff, max_diff);
    }
    return max_diff > m_max_height_difference;
  }

  void merge_segments(
    std::vector<Segment_2>& segments) {

    if (segments.size() == 0) return;

    const FT angle_threshold = FT(1);
    Segment_merger merger(
      m_ordinate_bound_2, angle_threshold);
    merger.merge_segments(segments);
  }

  void create_walls(
    const std::vector<Segment_2>& segments,
    std::vector<Approximate_face>& building_inner_walls) {

    if (segments.size() == 0) return;

    building_inner_walls.clear();
    building_inner_walls.reserve(segments.size());

    std::vector<Boundary> stub;
    const Building_walls_estimator westimator(
      stub, m_bottom_z, m_top_z);

    Approximate_face wall; Boundary boundary;
    for (const auto& segment : segments) {
      boundary.segment = segment;
      westimator.estimate_wall(boundary, wall.polygon);
      building_inner_walls.push_back(wall);
    }
  }
};

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_INTERNAL_PARTITION_23_ADAPTER_H
