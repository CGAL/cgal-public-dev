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

#ifndef CGAL_LEVELS_OF_DETAIL_INTERNAL_ESTIMATE_TREE_MODELS_H
#define CGAL_LEVELS_OF_DETAIL_INTERNAL_ESTIMATE_TREE_MODELS_H

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
  class Estimate_tree_models {

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
    Estimate_tree_models(
      const std::vector<Iterators>& clusters, 
      const Point_map& point_map) : 
    m_clusters(clusters), 
    m_point_map(point_map)
    { }

    void estimate(
      const FT min_radius_2, 
      std::vector<Tree_model>& models) const {

      std::vector<Tree_model> tmp_models;
      estimate_center_and_radius(tmp_models);
      clean_models(min_radius_2, tmp_models, models);
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

#endif // CGAL_LEVELS_OF_DETAIL_INTERNAL_ESTIMATE_TREE_MODELS_H
