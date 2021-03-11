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

#ifndef CGAL_LEVELS_OF_DETAIL_INTERNAL_SHAPE_DETECTION_LEAST_SQUARES_LINE_FIT_SORTING_H
#define CGAL_LEVELS_OF_DETAIL_INTERNAL_SHAPE_DETECTION_LEAST_SQUARES_LINE_FIT_SORTING_H

#include <CGAL/license/Levels_of_detail.h>

// STL includes.
#include <vector>
#include <algorithm>

// CGAL includes.
#include <CGAL/assertions.h>
#include <CGAL/Cartesian_converter.h>
#include <CGAL/linear_least_squares_fitting_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

// Internal includes.
#include <CGAL/Levels_of_detail/internal/utils.h>
#include <CGAL/Levels_of_detail/internal/property_map.h>

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

  template<
  typename GeomTraits,
  typename InputRange,
  typename NeighborQuery,
  typename PointMap>
  class Least_squares_line_fit_sorting {

  public:
    using Traits = GeomTraits;
    using Input_range = InputRange;
    using Neighbor_query = NeighborQuery;
    using Point_map = PointMap;
    using Seed_map = internal::Seed_property_map;

    Least_squares_line_fit_sorting(
      const Input_range& input_range,
      Neighbor_query& neighbor_query,
      const Point_map& point_map) :
    m_input_range(input_range),
    m_neighbor_query(neighbor_query),
    m_point_map(point_map) {

      CGAL_precondition(m_input_range.size() > 0);

      m_order.resize(m_input_range.size());
      for (std::size_t i = 0; i < m_input_range.size(); ++i)
        m_order[i] = i;
      m_scores.resize(m_input_range.size());
    }

    void sort() {

      compute_scores();
      CGAL_postcondition(m_scores.size() > 0);

      Compare_scores cmp(m_scores);
      std::sort(m_order.begin(), m_order.end(), cmp);
    }

    Seed_map seed_map() {
      return Seed_map(m_order);
    }

  private:
    using Local_traits = Exact_predicates_inexact_constructions_kernel;
    using Local_FT = typename Local_traits::FT;
    using Local_point_2 = typename Local_traits::Point_2;
    using Local_line_2 = typename Local_traits::Line_2;
    using To_local_converter = Cartesian_converter<Traits, Local_traits>;
    using Compare_scores = internal::Compare_scores<Local_FT>;

    void compute_scores() {

      std::vector<std::size_t> neighbors;
      std::vector<Local_point_2> points;

      for (std::size_t i = 0; i < m_input_range.size(); ++i) {

        neighbors.clear();
        m_neighbor_query(i, neighbors);
        neighbors.push_back(i);

        points.clear();
        for (std::size_t j = 0; j < neighbors.size(); ++j) {

          CGAL_precondition(neighbors[j] >= 0);
          CGAL_precondition(neighbors[j] < m_input_range.size());

          const auto& key = *(m_input_range.begin() + neighbors[j]);
          points.push_back(m_to_local_converter(get(m_point_map, key)));
        }
        CGAL_postcondition(points.size() == neighbors.size());

        Local_line_2  fitted_line;
        Local_point_2 fitted_centroid;

        m_scores[i] = linear_least_squares_fitting_2(
          points.begin(), points.end(),
          fitted_line, fitted_centroid,
          CGAL::Dimension_tag<0>());
      }
    }

    const Input_range& m_input_range;
    Neighbor_query& m_neighbor_query;
    const Point_map& m_point_map;

    std::vector<std::size_t> m_order;
    std::vector<Local_FT> m_scores;

    const To_local_converter m_to_local_converter;
  };

} // namespace internal
} // namespace Levels_of_detail
} // namespace CGAL

#endif // CGAL_LEVELS_OF_DETAIL_INTERNAL_SHAPE_DETECTION_LEAST_SQUARES_LINE_FIT_SORTING_H
