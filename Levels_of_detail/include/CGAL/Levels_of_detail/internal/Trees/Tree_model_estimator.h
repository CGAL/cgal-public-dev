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

#ifndef CGAL_LEVELS_OF_DETAIL_INTERNAL_TREE_MODEL_ESTIMATOR_H
#define CGAL_LEVELS_OF_DETAIL_INTERNAL_TREE_MODEL_ESTIMATOR_H

#include <CGAL/license/Levels_of_detail.h>

// STL includes.
#include <vector>
#include <utility>
#include <algorithm>

// Boost includes.
#include <boost/iterator/transform_iterator.hpp>

// CGAL includes.
#include <CGAL/barycenter.h>
#include <CGAL/number_utils.h>

// Internal includes.
#include <CGAL/Levels_of_detail/internal/struct.h>

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

  template<
  typename GeomTraits,
  typename InputRange,
  typename PointMap>
  class Tree_model_estimator {

  public:
    using Traits = GeomTraits;
    using Input_range = InputRange;
    using Point_map = PointMap;

    using FT = typename Traits::FT;
    using Point_2 = typename Traits::Point_2;
    using Point_3 = typename Traits::Point_3;
    using Circle_2 = typename Traits::Circle_2;

    using Iterator = typename Input_range::const_iterator;
    using Iterators = std::vector<Iterator>;
    using Tree_model = internal::Tree_model<Traits>;

  private:
    const std::vector<Iterators>& m_clusters;
    const Point_map& m_point_map;

    struct Point_from_iterator_and_pmap {
      using argument_type = Iterator;
      using return_type = std::pair<Point_3, FT>;

      const Point_map& m_point_map;
      Point_from_iterator_and_pmap(const Point_map& point_map) :
      m_point_map(point_map)
      { }

      return_type operator()(const argument_type& arg) const {
        return std::make_pair(get(m_point_map, *arg), FT(1));
      }
    };

  public:
    Tree_model_estimator(
      const std::vector<Iterators>& clusters,
      const Point_map& point_map) :
    m_clusters(clusters),
    m_point_map(point_map)
    { }

    void estimate_model_parameters(
      const FT min_radius_2,
      std::vector<Tree_model>& models) const {

      std::vector<Tree_model> tmp_models;
      estimate_center_and_radius(tmp_models);
      clean_models(min_radius_2, tmp_models, models);
    }

    void estimate_crown_parameters(
      std::vector<Tree_model>& models) const {

      for (auto& model : models) {
        const Point_2& center = model.center;
        const FT radius = model.radius;
        const std::size_t cluster_index = model.cluster_index;

        // Sort points.
        Iterators cluster = m_clusters[cluster_index];
        std::sort(cluster.begin(), cluster.end(),
        [&](const Iterator& a, const Iterator& b) -> bool {
          return get(m_point_map, *a).z() < get(m_point_map, *b).z();
        });

        // Set crown heights.
        const std::size_t size = cluster.size();
        const FT bn = get(m_point_map, *cluster[0]).z();
        const FT tn = get(m_point_map, *cluster[size-1]).z();
        const FT hn = tn - bn;
        const FT val = bn + hn / FT(4);

        std::size_t n = 0;
        for (std::size_t i = 0; i < cluster.size(); ++i) {
          const FT z = get(m_point_map, *cluster[i]).z();
          if (z > val) {
            n = i; break; }
        }

        const std::size_t m = size - n;
        const std::size_t idx0 = n;
        const std::size_t idx1 = std::size_t(n + m / 10.0);
        const std::size_t idx2 = std::size_t(n + m / 2.0);
        const std::size_t idx3 = std::size_t(n + 9.0 * m / 10.0);
        const std::size_t idx4 = size - 1;

        model.crown_z[0] = get(m_point_map, *cluster[idx0]).z();
        model.crown_z[1] = get(m_point_map, *cluster[idx1]).z();
        model.crown_z[2] = get(m_point_map, *cluster[idx2]).z();
        model.crown_z[3] = get(m_point_map, *cluster[idx3]).z();
        model.crown_z[4] = get(m_point_map, *cluster[idx4]).z();

        std::vector<FT> width{FT(0), FT(0), FT(0), FT(0)};
        std::vector<std::size_t> nb{0, 0, 0, 0};

        // Compute crown widths.
        for (const auto& it : cluster) {
          const Point_3& p = get(m_point_map, *it);
          const Point_2 q = Point_2(p.x(), p.y());

          std::size_t idx = 0;
          if (p.z() < model.crown_z[1])
            idx = 0;
          else if (p.z() < model.crown_z[2])
            idx = 1;
          else if (p.z() < model.crown_z[3])
            idx = 2;
          else idx = 3;

          width[idx] += CGAL::squared_distance(q, center);
          nb[idx]++;
        }

        // Normalize crown widths.
        for (std::size_t i = 0; i < width.size(); ++i)
          if (nb[i] != 0) width[i] = static_cast<FT>(
              CGAL::sqrt(CGAL::to_double(width[i] / static_cast<FT>(nb[i]))));

        // Set crown radiuses.
        model.crown_r[0] = model.trunk2_radius();
        model.crown_r[1] = (width[0] + width[1]) / FT(2);
        model.crown_r[2] = (width[1] + width[2]) / FT(2);
        model.crown_r[3] = (width[2] + width[3]) / FT(2);
      }
    }

  private:
    void estimate_center_and_radius(
      std::vector<Tree_model>& models) const {

      models.clear();
      models.reserve(m_clusters.size());
      std::size_t cluster_index = 0;
      for (const auto& cluster : m_clusters) {

        const Point_3 center = CGAL::barycenter(
          boost::make_transform_iterator(
            cluster.begin(), Point_from_iterator_and_pmap(m_point_map)),
          boost::make_transform_iterator(
            cluster.end(), Point_from_iterator_and_pmap(m_point_map)));

        FT radius = FT(0);
        for (const auto& it : cluster)
          radius += CGAL::squared_distance(center, get(m_point_map, *it));
        radius = static_cast<FT>(
          CGAL::sqrt(
            CGAL::to_double(
              radius / static_cast<FT>(cluster.size()))));

        models.push_back(
          Tree_model(Point_2(center.x(), center.y()), radius, cluster_index));
        ++cluster_index;
      }
      CGAL_assertion(models.size() == m_clusters.size());
    }

    void clean_models(
      const FT min_radius_2,
      std::vector<Tree_model>& input,
      std::vector<Tree_model>& output) const {

      output.clear();
      std::sort(input.begin(), input.end(),
        [](const Tree_model& a, const Tree_model& b) -> bool {
          return a.radius > b.radius;
        });

      for (const auto& a : input) {
        if (a.radius < min_radius_2)
          continue;

        bool okay = true;
        Circle_2 circle_a(a.center, a.radius * a.radius);
        for (const auto& b : output) {
          Circle_2 circle_b(b.center, b.radius * b.radius);

          if (CGAL::do_intersect(circle_a, circle_b) ||
            circle_b.has_on_bounded_side(circle_a.center())) {
            okay = false; break;
          }
        }
        if (okay)
          output.push_back(a);
      }
    }
  };

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_INTERNAL_TREE_MODEL_ESTIMATOR_H
