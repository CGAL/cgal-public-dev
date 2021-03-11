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

#ifndef CGAL_LEVELS_OF_DETAIL_INTERNAL_SHAPE_DETECTION_OVERLAPPING_REGION_H
#define CGAL_LEVELS_OF_DETAIL_INTERNAL_SHAPE_DETECTION_OVERLAPPING_REGION_H

#include <CGAL/license/Levels_of_detail.h>

// STL includes.
#include <vector>
#include <set>
#include <map>

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
  typename NormalMap,
  typename DefaultRegion>
  class Overlapping_region {

  public:
    using Traits = GeomTraits;
    using Input_range = InputRange;
    using Point_map = PointMap;
    using Normal_map = NormalMap;
    using Region = DefaultRegion;

    using FT = typename Traits::FT;
    using Point_3 = typename Traits::Point_3;
    using Vector_3 = typename Traits::Vector_3;
    using Plane_3 = typename Traits::Plane_3;

    using Indices = std::vector<std::size_t>;
    using Groups = std::vector< std::set<std::size_t> >;

    using Squared_length_3 = typename Traits::Compute_squared_length_3;
    using Squared_distance_3 = typename Traits::Compute_squared_distance_3;
    using Scalar_product_3 = typename Traits::Compute_scalar_product_3;

    using Get_sqrt = internal::Get_sqrt<Traits>;
    using Sqrt = typename Get_sqrt::Sqrt;

    using Size_pair = std::pair<std::size_t, std::size_t>;

    Overlapping_region(
      const Input_range& input_range,
      const std::vector<Indices>& roofs,
      Region& region,
      const FT distance_threshold,
      const FT angle_threshold,
      const Point_map& point_map,
      const Normal_map& normal_map,
      const Traits traits = Traits()) :
    m_input_range(input_range),
    m_roofs(roofs),
    m_region(region),
    m_distance_threshold(distance_threshold),
    m_normal_threshold(static_cast<FT>(
      std::cos(
        CGAL::to_double(
          (angle_threshold * static_cast<FT>(CGAL_PI)) / FT(180))))),
    m_point_map(point_map),
    m_normal_map(normal_map),
    m_squared_length_3(traits.compute_squared_length_3_object()),
    m_squared_distance_3(traits.compute_squared_distance_3_object()),
    m_scalar_product_3(traits.compute_scalar_product_3_object()),
    m_sqrt(Get_sqrt::sqrt_object(traits)) {

      CGAL_precondition(m_input_range.size() > 0);

      CGAL_precondition(m_distance_threshold >= FT(0));
      CGAL_precondition(m_normal_threshold >= FT(0) && m_normal_threshold <= FT(1));

      m_roof_map.clear();
      for (std::size_t i = 0; i < m_input_range.size(); ++i)
        m_roof_map[i] = std::size_t(-1);
      for (std::size_t k = 0; k < m_roofs.size(); ++k)
        for (const std::size_t idx : m_roofs[k])
          m_roof_map[idx] = k;

      Input_range range;
      for (std::size_t k = 0; k < m_roofs.size(); ++k) {
        const auto& roof = m_roofs[k];

        range.clear();
        for (const std::size_t idx : roof)
          range.push_back(m_input_range[idx]);

        Plane_3 plane;
        internal::plane_from_points_3(range, m_point_map, plane);
        m_planes[k] = plane;

        Vector_3 normal = plane.orthogonal_vector();
        const FT normal_length = m_sqrt(m_squared_length_3(normal));
        normal /= normal_length;
        m_normals[k] = normal;
      }

      m_groups.clear();
      m_groups.resize(m_input_range.size());
      for (std::size_t i = 0; i < m_input_range.size(); ++i) {
        const std::size_t roof_idx = m_roof_map.at(i);
        if (roof_idx != std::size_t(-1))
          m_groups[i].insert(roof_idx);
      }
    }

    bool is_already_visited(
      const std::size_t,
      const std::size_t query_index,
      const bool is_visited) const {

      const std::size_t roof_idx = m_roof_map.at(query_index);
      if (
        !is_visited &&
        roof_idx != std::size_t(-1) &&
        m_current_roof_idx != std::size_t(-1) &&
        roof_idx == m_current_roof_idx)
        return true;
      return false;
    }

    bool is_part_of_region(
      const std::size_t stub,
      const std::size_t query_index,
      const std::vector<std::size_t>& region) {

      if (m_current_roof_idx == std::size_t(-1))
        return false;

      const std::size_t roof_idx = m_roof_map.at(query_index);
      if (roof_idx == std::size_t(-1))
        return false;

      if (m_current_roof_idx == roof_idx)
        return false;

      const auto& plane_of_best_fit  = m_planes.at(m_current_roof_idx);
      const auto& normal_of_best_fit = m_normals.at(m_current_roof_idx);

      const auto& key = *(m_input_range.begin() + query_index);
      const Point_3& query_point = get(m_point_map, key);

      const Vector_3& normal = get(m_normal_map, key);
      const FT normal_length = m_sqrt(m_squared_length_3(normal));
      const Vector_3 query_normal = normal / normal_length;

      const FT distance_to_fitted_plane =
      m_sqrt(m_squared_distance_3(query_point, plane_of_best_fit));

      const FT cos_value =
      CGAL::abs(m_scalar_product_3(query_normal, normal_of_best_fit));

      const bool result = (( distance_to_fitted_plane <= m_distance_threshold ) &&
        ( cos_value >= m_normal_threshold ));

      if (result)
        m_groups[query_index].insert(m_current_roof_idx);
      return result;
    }

    bool is_valid_region(const std::vector<std::size_t>& region) const {
      return m_region.is_valid_region(region);
    }

    void update(const std::vector<std::size_t>& region) {

      CGAL_precondition(region.size() > 0);
      if (region.size() == 1)
        m_current_roof_idx = m_roof_map.at(region[0]);
    }

    void get_overlapping_points(
      std::map<std::size_t, bool>& overlapping) {

      overlapping.clear();
      std::size_t count = 0;
      for (std::size_t i = 0; i < m_groups.size(); ++i) {
        if (m_groups[i].size() > 1) {
          ++count; overlapping[i] = true;
        } else
          overlapping[i] = false;
      }
      std::cout << "Num overlapping points: " << count << std::endl;
    }

    void get_groups(Groups& groups) {
      groups = m_groups;
    }

  private:
    const Input_range& m_input_range;
    const std::vector<Indices>& m_roofs;

    Region& m_region;

    const FT m_distance_threshold;
    const FT m_normal_threshold;

    const Point_map& m_point_map;
    const Normal_map& m_normal_map;

    const Squared_length_3 m_squared_length_3;
    const Squared_distance_3 m_squared_distance_3;
    const Scalar_product_3 m_scalar_product_3;
    const Sqrt m_sqrt;

    std::map<std::size_t, std::size_t> m_roof_map;

    std::map<std::size_t, Plane_3> m_planes;
    std::map<std::size_t, Vector_3> m_normals;

    std::size_t m_current_roof_idx;

    Groups m_groups;
  };

} // namespace internal
} // namespace Levels_of_detail
} // namespace CGAL

#endif // CGAL_LEVELS_OF_DETAIL_INTERNAL_SHAPE_DETECTION_OVERLAPPING_REGION_H
