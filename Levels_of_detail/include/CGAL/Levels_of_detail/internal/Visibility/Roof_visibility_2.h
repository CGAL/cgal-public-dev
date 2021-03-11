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

#ifndef CGAL_LEVELS_OF_DETAIL_INTERNAL_ROOF_VISIBILITY_2_H
#define CGAL_LEVELS_OF_DETAIL_INTERNAL_ROOF_VISIBILITY_2_H

// STL includes.
#include <set>
#include <vector>
#include <utility>
#include <iostream>
#include <fstream>
#include <memory>
#include <algorithm>

// CGAL includes.
#include <CGAL/assertions.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/Random.h>

// Internal includes.
#include <CGAL/Levels_of_detail/enum.h>
#include <CGAL/Levels_of_detail/internal/struct.h>
#include <CGAL/Levels_of_detail/internal/utils.h>

// Spatial search.
#include <CGAL/Levels_of_detail/internal/Spatial_search/K_neighbor_query.h>

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

  template<typename GeomTraits>
  class Roof_visibility_2 {

  public:
    using Traits = GeomTraits;

    using FT = typename Traits::FT;
    using Point_2 = typename Traits::Point_2;
    using Point_3 = typename Traits::Point_3;
    using Triangle_2 = typename Traits::Triangle_2;

    using Indices = std::vector<std::size_t>;
    using Partition_2 = internal::Partition_2<Traits>;
    using Face = typename Partition_2::Face;
    using Building = internal::Building<Traits>;

    using Point_generator = CGAL::Random_points_in_triangle_2<Point_2>;
    using Pair = std::pair<Point_2, std::size_t>;
    using Point_map_2 = CGAL::First_of_pair_property_map<Pair>;
    using Point_map_3 = CGAL::Identity_property_map<Point_3>;
    using K_neighbor_query = internal::K_neighbor_query<Traits, std::vector<Pair>, Point_map_2>;
    using Random = CGAL::Random;
    using Triangulation = Triangulation<Traits>;
    using Location_type = typename Triangulation::Delaunay::Locate_type;

    Roof_visibility_2(
      const std::vector<Point_3>& input_range,
      const Building& building,
      const std::vector<Indices>& roof_points_3) :
    m_input_range(input_range),
    m_building(building),
    m_roof_points_3(roof_points_3),
    m_num_samples(100), // num samples per triangle
    m_k(1),
    m_random(0)
    { }

    void compute(Partition_2& partition) {

      compute_in_out_visibility(partition);
      compute_roof_visibility(partition);
    }

    const std::size_t number_of_actual_roofs(
      const Partition_2& partition) const {
      return compute_number_of_actual_roofs(partition);
    }

  private:
    const std::vector<Point_3>& m_input_range;
    const Building& m_building;
    const std::vector<Indices>& m_roof_points_3;

    const std::size_t m_num_samples;
    std::vector<Point_2> m_samples;
    const std::size_t m_k;
    Random m_random;

    Point_map_2 m_point_map_2;
    Point_map_3 m_point_map_3;

    std::vector<std::size_t> m_roof_indices;
    std::vector<Pair> m_queries;
    std::shared_ptr<K_neighbor_query> m_neighbor_query_ptr;

    std::size_t compute_number_of_actual_roofs(
      const Partition_2& partition) const {

      std::set<std::size_t> roofs;
      for (const auto& pface : partition.faces) {
        if (pface.visibility == Visibility_label::OUTSIDE)
          continue;
        if (pface.label == std::size_t(-1))
          continue;
        roofs.insert(pface.label);
      }
      return roofs.size();
    }

    void compute_in_out_visibility(Partition_2& partition) {
      const auto& ref = m_building.base1.triangulation.delaunay;

      for (auto& face : partition.faces) {
        const auto& polygon = face.outer_polygon;

        if (polygon.size() < 3) {

          face.visibility = Visibility_label::OUTSIDE;
          face.inside = FT(0); face.outside = FT(1); continue;
        }

        const auto& p0 = polygon[0]; std::size_t count = 0;
        for (std::size_t i = 1; i < polygon.size() - 1; ++i) {
          const std::size_t ip = i + 1;

          const auto& p1 = polygon[i];
          const auto& p2 = polygon[ip];

          const FT x = (p0.x() + p1.x() + p2.x()) / FT(3);
          const FT y = (p0.y() + p1.y() + p2.y()) / FT(3);

          const Point_2 b = Point_2(x, y);
          Location_type type; int stub;
          const auto bh = ref.locate(b, type, stub);

          const Triangle_2 triangle = Triangle_2(p0, p1, p2);
          const FT area = CGAL::abs(triangle.area());

          if (
            bh->info().tagged &&
            type == Triangulation::Delaunay::FACE &&
            area > FT(1) / FT(1000))
            ++count;
        }

        if (count >= 1) {
          face.visibility = Visibility_label::INSIDE;
          face.inside = FT(1); face.outside = FT(0);
        } else {
          face.visibility = Visibility_label::OUTSIDE;
          face.inside = FT(0); face.outside = FT(1);
        }
      }
    }

    void compute_roof_visibility(Partition_2& partition) {

      create_tree();
      for (auto& face : partition.faces)
        if (face.visibility == Visibility_label::INSIDE)
          compute_face_label(face);
    }

    void create_tree() {

      std::size_t num_points = 0;
      for (const auto& region : m_roof_points_3)
        num_points += region.size();

      m_queries.clear();
      m_queries.reserve(num_points);

      m_roof_indices.clear();
      m_roof_indices.reserve(num_points);

      for (std::size_t i = 0; i < m_roof_points_3.size(); ++i) {
        for (const std::size_t idx : m_roof_points_3[i]) {
          const Point_3& p = get(m_point_map_3, *(m_input_range.begin() + idx));
          const Point_2 q = internal::point_2_from_point_3(p);
          m_queries.push_back(std::make_pair(q, i));
          m_roof_indices.push_back(i);
        }
      }

      m_neighbor_query_ptr = std::make_shared<K_neighbor_query>(
        m_queries, FT(m_k), m_point_map_2);
    }

    void compute_face_label(Face& face) {

      Indices neighbors;
      const auto& polygon = face.outer_polygon;

      std::vector<int> max_count(
        m_roof_points_3.size(), 0);

      const auto& p0 = polygon[0];
      for (std::size_t i = 1; i < polygon.size() - 1; ++i) {
        const std::size_t ip = i + 1;

        const auto& p1 = polygon[i];
        const auto& p2 = polygon[ip];
        const Triangle_2 triangle = Triangle_2(p0, p1, p2);

        const FT area = CGAL::abs(triangle.area());
        if (area <= FT(1) / FT(1000))
          continue;

        m_samples.clear();
        Point_generator generator(triangle, m_random);
        std::copy_n(
          generator, m_num_samples, std::back_inserter(m_samples));

        for (const auto& p : m_samples) {
          (*m_neighbor_query_ptr)(p, neighbors);

          for (const std::size_t idx : neighbors) {
            const std::size_t label = m_queries[idx].second;
            max_count[label] += 1;
          }
        }
      }

      std::size_t face_label = std::size_t(-1);
      int max_value = -1;
      for (std::size_t i = 0; i < max_count.size(); ++i) {
        const int value = max_count[i];
        if (value == 0) continue;

        if (value > max_value) {
          max_value = value; face_label = i;
        }
      }

      face.label = face_label;
      face.probabilities.clear();
      face.probabilities.resize(max_count.size(), FT(0));

      if (max_value != -1) {
        for (std::size_t i = 0; i < max_count.size(); ++i) {

          const FT probability =
            static_cast<FT>(max_count[i]) / static_cast<FT>(max_value);
          face.probabilities[i] = probability;
        }
      }
    }
  };

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_INTERNAL_ROOF_VISIBILITY_2_H
