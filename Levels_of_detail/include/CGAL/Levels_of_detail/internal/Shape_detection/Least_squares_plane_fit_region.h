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

#ifndef CGAL_LEVELS_OF_DETAIL_INTERNAL_SHAPE_DETECTION_LEAST_SQUARES_PLANE_FIT_REGION_H
#define CGAL_LEVELS_OF_DETAIL_INTERNAL_SHAPE_DETECTION_LEAST_SQUARES_PLANE_FIT_REGION_H

#include <CGAL/license/Levels_of_detail.h>

// STL includes.
#include <vector>

// CGAL includes.
#include <CGAL/assertions.h>
#include <CGAL/number_utils.h>
#include <CGAL/Cartesian_converter.h>
#include <CGAL/linear_least_squares_fitting_3.h>
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
  class Least_squares_plane_fit_region {

  public:
    using Traits = GeomTraits;
    using Input_range = InputRange;
    using Point_map = PointMap;
    using Normal_map = NormalMap;

    using FT = typename Traits::FT;
    using Point_3 = typename Traits::Point_3;
    using Vector_3 = typename Traits::Vector_3;
    using Plane_3 = typename Traits::Plane_3;

    using Local_traits = Exact_predicates_inexact_constructions_kernel;
    using Local_point_3 = typename Local_traits::Point_3;
    using Local_plane_3 = typename Local_traits::Plane_3;
    using To_local_converter = Cartesian_converter<Traits, Local_traits>;

    using Squared_length_3 = typename Traits::Compute_squared_length_3;
    using Squared_distance_3 = typename Traits::Compute_squared_distance_3;
    using Scalar_product_3 = typename Traits::Compute_scalar_product_3;

    using Get_sqrt = internal::Get_sqrt<Traits>;
    using Sqrt = typename Get_sqrt::Sqrt;

    Least_squares_plane_fit_region(
      const Input_range& input_range,
      const FT distance_threshold,
      const FT angle_threshold,
      const FT min_area,
      const FT distance_to_line,
      const FT alpha,
      const Point_map& point_map,
      const Normal_map& normal_map,
      const Traits traits = Traits()) :
    m_input_range(input_range),
    m_distance_threshold(distance_threshold),
    m_normal_threshold(static_cast<FT>(
      std::cos(
        CGAL::to_double(
          (angle_threshold * static_cast<FT>(CGAL_PI)) / FT(180))))),
    m_min_area(min_area),
    m_distance_to_line(distance_to_line),
    m_alpha(alpha),
    m_point_map(point_map),
    m_normal_map(normal_map),
    m_squared_length_3(traits.compute_squared_length_3_object()),
    m_squared_distance_3(traits.compute_squared_distance_3_object()),
    m_scalar_product_3(traits.compute_scalar_product_3_object()),
    m_sqrt(Get_sqrt::sqrt_object(traits)) {

      CGAL_precondition(m_input_range.size() > 0);

      CGAL_precondition(m_distance_threshold >= FT(0));
      CGAL_precondition(m_normal_threshold >= FT(0) && m_normal_threshold <= FT(1));
      CGAL_precondition(m_min_area > FT(0));
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
      const Point_3& query_point = get(m_point_map, key);

      const Vector_3& normal = get(m_normal_map, key);
      const FT normal_length = m_sqrt(m_squared_length_3(normal));
      CGAL_precondition(normal_length > FT(0));
      const Vector_3 query_normal = normal / normal_length;

      const FT distance_to_fitted_plane =
      m_sqrt(m_squared_distance_3(query_point, m_plane_of_best_fit));

      const FT cos_value =
      CGAL::abs(m_scalar_product_3(query_normal, m_normal_of_best_fit));

      return (( distance_to_fitted_plane <= m_distance_threshold ) &&
        ( cos_value >= m_normal_threshold ));
    }

    bool is_valid_region(const std::vector<std::size_t>& region) const {

      if (region.size() < 3) return false;
      const FT distance = internal::average_distance_to_line_3(
      m_input_range, m_point_map, region);
      const FT area = internal::points_area(
      m_input_range, m_point_map, region, m_alpha);
      return ( distance >= m_distance_to_line && area >= m_min_area );
    }

    void update(const std::vector<std::size_t>& region) {

      CGAL_precondition(region.size() > 0);
      if (region.size() == 1) { // create new reference plane and normal

        CGAL_precondition(region[0] >= 0);
        CGAL_precondition(region[0] < m_input_range.size());

        // The best fit plane will be a plane through this point with
        // its normal being the point's normal.
        const auto& key = *(m_input_range.begin() + region[0]);

        const Point_3& point = get(m_point_map, key);
        const Vector_3& normal = get(m_normal_map, key);

        const FT normal_length = m_sqrt(m_squared_length_3(normal));

        CGAL_precondition(normal_length > FT(0));
        m_normal_of_best_fit =
        normal / normal_length;

        m_plane_of_best_fit =
        Plane_3(point, m_normal_of_best_fit);

      } else { // update reference plane and normal

        std::vector<Local_point_3> points;
        points.reserve(region.size());

        for (std::size_t i = 0; i < region.size(); ++i) {

          CGAL_precondition(region[i] >= 0);
          CGAL_precondition(region[i] < m_input_range.size());

          const auto& key = *(m_input_range.begin() + region[i]);
          points.push_back(m_to_local_converter(get(m_point_map, key)));
        }
        CGAL_postcondition(points.size() == region.size());

        Local_plane_3 fitted_plane;
        Local_point_3 fitted_centroid;

        // The best fit plane will be a plane fitted to all region points with
        // its normal being perpendicular to the plane.
        CGAL::linear_least_squares_fitting_3(
          points.begin(), points.end(),
          fitted_plane, fitted_centroid,
          CGAL::Dimension_tag<0>());

        m_plane_of_best_fit =
        Plane_3(
          static_cast<FT>(fitted_plane.a()),
          static_cast<FT>(fitted_plane.b()),
          static_cast<FT>(fitted_plane.c()),
          static_cast<FT>(fitted_plane.d()));

        const Vector_3 normal = m_plane_of_best_fit.orthogonal_vector();
        const FT normal_length = m_sqrt(m_squared_length_3(normal));

        CGAL_precondition(normal_length > FT(0));
        m_normal_of_best_fit = normal / normal_length;
      }
    }

  private:
    const Input_range& m_input_range;

    const FT m_distance_threshold;
    const FT m_normal_threshold;
    const FT m_min_area;
    const FT m_distance_to_line;
    const FT m_alpha;

    const Point_map& m_point_map;
    const Normal_map& m_normal_map;

    const Squared_length_3 m_squared_length_3;
    const Squared_distance_3 m_squared_distance_3;
    const Scalar_product_3 m_scalar_product_3;
    const Sqrt m_sqrt;

    const To_local_converter m_to_local_converter;

    Plane_3 m_plane_of_best_fit;
    Vector_3 m_normal_of_best_fit;
  };

} // namespace internal
} // namespace Levels_of_detail
} // namespace CGAL

#endif // CGAL_LEVELS_OF_DETAIL_INTERNAL_SHAPE_DETECTION_LEAST_SQUARES_PLANE_FIT_REGION_H
