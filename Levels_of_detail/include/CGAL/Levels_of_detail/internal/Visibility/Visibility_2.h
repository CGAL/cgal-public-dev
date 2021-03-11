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

#ifndef CGAL_LEVELS_OF_DETAIL_INTERNAL_VISIBILITY_2_H
#define CGAL_LEVELS_OF_DETAIL_INTERNAL_VISIBILITY_2_H

// STL includes.
#include <vector>
#include <utility>

// CGAL includes.
#include <CGAL/assertions.h>
#include <CGAL/property_map.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/Random.h>

// LOD includes.
#include <CGAL/Levels_of_detail/enum.h>

// Internal includes.
#include <CGAL/Levels_of_detail/internal/utils.h>
#include <CGAL/Levels_of_detail/internal/struct.h>

// Spatial search.
#include <CGAL/Levels_of_detail/internal/Spatial_search/K_neighbor_query.h>
#include <CGAL/Levels_of_detail/internal/Spatial_search/Sphere_neighbor_query.h>

// Testing.
#include "../../../../../test/Levels_of_detail/include/Saver.h"

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

template<
typename GeomTraits,
typename InputRange,
typename PointMap2,
typename VisibilityMap>
class Visibility_2 {

public:
  using Traits = GeomTraits;
  using Input_range = InputRange;
  using Point_map_2 = PointMap2;
  using Visibility_map = VisibilityMap;

  using FT = typename Traits::FT;
  using Point_2 = typename Traits::Point_2;
  using Point_3 = typename Traits::Point_3;
  using Triangle_2 = typename Traits::Triangle_2;

  using Points_2 = std::vector<Point_2>;
  using Indices = std::vector<std::size_t>;

  using Triangulation = typename internal::Triangulation<Traits>::Delaunay;
  using Location_type = typename Triangulation::Locate_type;

  using Partition_2 = internal::Partition_2<Traits>;

  using Saver = Saver<Traits>;
  using Color = CGAL::Color;

  using Random = CGAL::Random;
  using Point_generator = CGAL::Random_points_in_triangle_2<Point_2>;

  using Identity_map_2 = CGAL::Identity_property_map<Point_2>;

  using K_neighbor_query =
    internal::K_neighbor_query<Traits, Input_range, Point_map_2>;
  using Sphere_neighbor_query =
    internal::Sphere_neighbor_query<Traits, Input_range, Point_map_2>;

  Visibility_2(
    const Input_range& input_range,
    const Point_map_2& point_map_2,
    const Visibility_map& visibility_map,
    const FT threshold = FT(1) / FT(2)) :
  m_input_range(input_range),
  m_point_map_2(point_map_2),
  m_visibility_map(visibility_map),
  m_random(0),
  m_k(FT(1)), // for k neighbor query
  m_radius(FT(1) / FT(4)), // for sphere neighbor query
  m_samples_per_triangle(100),
  m_threshold(threshold)
  { }

  void compute(Partition_2& partition_2) {

    // Create tree.

    // Also change below!
    // using Neighbor_query = K_neighbor_query;
    // Neighbor_query neighbor_query(
    //   m_input_range, m_k, m_point_map_2);

    using Neighbor_query = Sphere_neighbor_query;
    Neighbor_query neighbor_query(
      m_input_range, m_radius, m_point_map_2);

    // Compute visibility.
    Indices neighbors;
    std::vector<Point_2> samples;

    for (auto& face : partition_2.faces) {
      const auto& polygon = face.outer_polygon;

      if (polygon.size() < 3) {

        face.visibility = Visibility_label::OUTSIDE;
        face.inside = FT(0); face.outside = FT(1); continue;
      }

      FT in = FT(1); FT out = FT(1);
      const auto& p0 = polygon[0];
      for (std::size_t i = 1; i < polygon.size() - 1; ++i) {
        const std::size_t ip = i + 1;

        const auto& p1 = polygon[i];
        const auto& p2 = polygon[ip];
        const Triangle_2 triangle = Triangle_2(p0, p1, p2);

        const FT area = CGAL::abs(triangle.area());
        if (area <= FT(1) / FT(1000))
          continue;

        samples.clear();
        Point_generator generator(triangle, m_random);
        std::copy_n(
          generator, m_samples_per_triangle, std::back_inserter(samples));

        for (const auto& p : samples) {
          neighbor_query(p, neighbors);

          // Use with sphere neighbor query!
          if (neighbors.size() == 0) out += FT(1);
          else in += FT(1);

          /* Use with k neighbor query!
          for (const std::size_t idx : neighbors) {
            const bool is_inside = get(m_visibility_map, *(m_input_range.begin() + idx));
            if (is_inside) in += FT(1);
            else out += FT(1);
          } */
        }
      }

      const FT sum = in + out;
      in /= sum; out /= sum;

      CGAL_assertion(in >= FT(0) && in <= FT(1));
      CGAL_assertion(out >= FT(0) && out <= FT(1));

      if (in > m_threshold) face.visibility = Visibility_label::INSIDE;
      else face.visibility = Visibility_label::OUTSIDE;
      face.inside = in; face.outside = out;
    }

    // save_samples(samples, "/Users/monet/Documents/gf/lod/logs/buildings/tmp/samples");
  }

private:
  const Input_range& m_input_range;
  const Point_map_2& m_point_map_2;
  const Visibility_map& m_visibility_map;

  Random m_random;

  const FT m_k;
  const FT m_radius;

  const std::size_t m_samples_per_triangle;
  const FT m_threshold;

  Saver m_saver;

  void save_samples(
    const std::vector<Point_2>& samples,
    const std::string name) {

    std::vector<Point_3> points;
    points.reserve(samples.size());
    for (const auto& p: samples)
      points.push_back(Point_3(p.x(), p.y(), FT(0)));
    const Color color(0, 0, 0);
    m_saver.export_points(points, color, name);
  }
};

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_INTERNAL_VISIBILITY_2_H
