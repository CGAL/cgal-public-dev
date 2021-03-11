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

#ifndef CGAL_LEVELS_OF_DETAIL_INTERNAL_SHAPE_DETECTION_LEAST_SQUARES_LINE_FIT_REGION_H
#define CGAL_LEVELS_OF_DETAIL_INTERNAL_SHAPE_DETECTION_LEAST_SQUARES_LINE_FIT_REGION_H

#include <CGAL/license/Levels_of_detail.h>

// STL includes.
#include <vector>

// CGAL includes.
#include <CGAL/assertions.h>
#include <CGAL/number_utils.h>
#include <CGAL/Cartesian_converter.h>
#include <CGAL/linear_least_squares_fitting_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

// Internal includes.
#include <CGAL/Levels_of_detail/internal/utils.h>

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

  template<
  typename GeomTraits,
  typename InputRange,
  typename PointMap,
  typename NormalMap>
  class Least_squares_line_fit_region {

  public:
    using Traits = GeomTraits;
    using Input_range = InputRange;
    using Point_map = PointMap;
    using Normal_map = NormalMap;

    using FT = typename Traits::FT;
    using Point_2 = typename Traits::Point_2;
    using Vector_2 = typename Traits::Vector_2;
    using Line_2 = typename Traits::Line_2;

    using Local_traits = Exact_predicates_inexact_constructions_kernel;
    using Local_point_2 = typename Local_traits::Point_2;
    using Local_line_2 = typename Local_traits::Line_2;
    using To_local_converter = Cartesian_converter<Traits, Local_traits>;

    using Squared_length_2 = typename Traits::Compute_squared_length_2;
    using Squared_distance_2 = typename Traits::Compute_squared_distance_2;
    using Scalar_product_2 = typename Traits::Compute_scalar_product_2;

    using Get_sqrt = internal::Get_sqrt<Traits>;
    using Sqrt = typename Get_sqrt::Sqrt;

    Least_squares_line_fit_region(
      const Input_range& input_range,
      const FT distance_threshold,
      const FT angle_threshold,
      const FT min_length,
      const Point_map& point_map,
      const Normal_map& normal_map,
      const Traits traits = Traits()) :
    m_input_range(input_range),
    m_distance_threshold(distance_threshold),
    m_normal_threshold(static_cast<FT>(
      std::cos(
        CGAL::to_double(
          (angle_threshold * static_cast<FT>(CGAL_PI)) / FT(180))))),
    m_min_length(min_length),
    m_point_map(point_map),
    m_normal_map(normal_map),
    m_squared_length_2(traits.compute_squared_length_2_object()),
    m_squared_distance_2(traits.compute_squared_distance_2_object()),
    m_scalar_product_2(traits.compute_scalar_product_2_object()),
    m_sqrt(Get_sqrt::sqrt_object(traits)) {

      CGAL_precondition(m_input_range.size() > 0);

      CGAL_precondition(m_distance_threshold >= FT(0));
      CGAL_precondition(m_normal_threshold >= FT(0) && m_normal_threshold <= FT(1));
      CGAL_precondition(m_min_length > FT(0));
    }

    bool is_already_visited(
      const std::size_t,
      const std::size_t query_index,
      const bool is_visited) const { return false; }

    bool is_part_of_region(
      const std::size_t,
      const std::size_t query_index,
      const std::vector<std::size_t>&) const {

      CGAL_precondition(query_index >= 0);
      CGAL_precondition(query_index < m_input_range.size());

      const auto& key = *(m_input_range.begin() + query_index);
      const Point_2& query_point = get(m_point_map, key);

      const Vector_2& normal = get(m_normal_map, key);
      const FT normal_length = m_sqrt(m_squared_length_2(normal));
      CGAL_precondition(normal_length > FT(0));
      const Vector_2 query_normal = normal / normal_length;

      const FT distance_to_fitted_line =
      m_sqrt(m_squared_distance_2(query_point, m_line_of_best_fit));

      const FT cos_value =
      CGAL::abs(m_scalar_product_2(query_normal, m_normal_of_best_fit));

      return (( distance_to_fitted_line <= m_distance_threshold ) &&
        ( cos_value >= m_normal_threshold ));
    }

    bool is_valid_region(const std::vector<std::size_t>& region) const {

      if (region.size() < 2) return false;
      const FT squared_min_length = m_min_length * m_min_length;
      return ( internal::points_squared_length_2(
        m_input_range, m_point_map, region, m_line_of_best_fit) >= squared_min_length );
    }

    void update(const std::vector<std::size_t>& region) {

      CGAL_precondition(region.size() > 0);
      if (region.size() == 1) { // create new reference line and normal

        CGAL_precondition(region[0] >= 0);
        CGAL_precondition(region[0] < m_input_range.size());

        // The best fit line will be a line through this point with
        // its normal being the point's normal.
        const auto& key = *(m_input_range.begin() + region[0]);

        const Point_2& point = get(m_point_map, key);
        const Vector_2& normal = get(m_normal_map, key);

        const FT normal_length = m_sqrt(m_squared_length_2(normal));

        CGAL_precondition(normal_length > FT(0));
        m_normal_of_best_fit =
        normal / normal_length;

        m_line_of_best_fit =
        Line_2(point, m_normal_of_best_fit).perpendicular(point);

      } else { // update reference line and normal

        std::vector<Local_point_2> points;
        points.reserve(region.size());

        for (std::size_t i = 0; i < region.size(); ++i) {

          CGAL_precondition(region[i] >= 0);
          CGAL_precondition(region[i] < m_input_range.size());

          const auto& key = *(m_input_range.begin() + region[i]);
          points.push_back(m_to_local_converter(get(m_point_map, key)));
        }
        CGAL_postcondition(points.size() == region.size());

        Local_line_2 fitted_line;
        Local_point_2 fitted_centroid;

        // The best fit line will be a line fitted to all region points with
        // its normal being perpendicular to the line.
        CGAL::linear_least_squares_fitting_2(
          points.begin(), points.end(),
          fitted_line, fitted_centroid,
          CGAL::Dimension_tag<0>());

        m_line_of_best_fit =
        Line_2(
          static_cast<FT>(fitted_line.a()),
          static_cast<FT>(fitted_line.b()),
          static_cast<FT>(fitted_line.c()));

        const Vector_2 normal =
        m_line_of_best_fit.perpendicular(m_line_of_best_fit.point(0)).to_vector();
        const FT normal_length = m_sqrt(m_squared_length_2(normal));

        CGAL_precondition(normal_length > FT(0));
        m_normal_of_best_fit = normal / normal_length;
      }
    }

  private:
    const Input_range& m_input_range;

    const FT m_distance_threshold;
    const FT m_normal_threshold;
    const FT m_min_length;

    const Point_map& m_point_map;
    const Normal_map& m_normal_map;

    const Squared_length_2 m_squared_length_2;
    const Squared_distance_2 m_squared_distance_2;
    const Scalar_product_2 m_scalar_product_2;
    const Sqrt m_sqrt;

    const To_local_converter m_to_local_converter;

    Line_2 m_line_of_best_fit;
    Vector_2 m_normal_of_best_fit;
  };

} // namespace internal
} // namespace Levels_of_detail
} // namespace CGAL

#endif // CGAL_LEVELS_OF_DETAIL_INTERNAL_SHAPE_DETECTION_LEAST_SQUARES_LINE_FIT_REGION_H
