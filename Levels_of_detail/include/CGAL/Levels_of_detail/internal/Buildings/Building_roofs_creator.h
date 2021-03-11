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

#ifndef CGAL_LEVELS_OF_DETAIL_INTERNAL_ROOFS_CREATOR_H
#define CGAL_LEVELS_OF_DETAIL_INTERNAL_ROOFS_CREATOR_H

// STL includes.
#include <vector>

// CGAL includes.
#include <CGAL/assertions.h>

// Internal includes.
#include <CGAL/Levels_of_detail/internal/utils.h>
#include <CGAL/Levels_of_detail/internal/struct.h>

// Spatial search.
#include <CGAL/Levels_of_detail/internal/Spatial_search/Sphere_neighbor_query.h>

// Shape detection.
#include <CGAL/Levels_of_detail/internal/Shape_detection/Region_growing.h>
#include <CGAL/Levels_of_detail/internal/Shape_detection/Estimate_normals_3.h>
#include <CGAL/Levels_of_detail/internal/Shape_detection/Least_squares_plane_fit_region.h>
#include <CGAL/Levels_of_detail/internal/Shape_detection/Least_squares_plane_fit_sorting.h>
#include <CGAL/Levels_of_detail/internal/Shape_detection/Overlapping_region.h>

// Testing.
#include "../../../../../test/Levels_of_detail/include/Saver.h"

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

  template<
  typename GeomTraits,
  typename PointMap3>
  class Building_roofs_creator {

  public:
    using Traits = GeomTraits;
    using Point_map_3 = PointMap3;

    using FT = typename Traits::FT;
    using Point_3 = typename Traits::Point_3;
    using Vector_3 = typename Traits::Vector_3;
    using Plane_3 = typename Traits::Plane_3;
    using Indices = std::vector<std::size_t>;

    using Pair_item_3 = std::pair<Point_3, Vector_3>;
    using Pair_range_3 = std::vector<Pair_item_3>;
    using First_of_pair_map = CGAL::First_of_pair_property_map<Pair_item_3>;
    using Second_of_pair_map = CGAL::Second_of_pair_property_map<Pair_item_3>;

    using Sphere_neighbor_query =
    internal::Sphere_neighbor_query<Traits, Indices, Point_map_3>;
    using Normal_estimator_3 =
    internal::Estimate_normals_3<Traits, Indices, Point_map_3, Sphere_neighbor_query>;
    using LSPF_region =
    internal::Least_squares_plane_fit_region<Traits, Pair_range_3, First_of_pair_map, Second_of_pair_map>;
    using OLPF_region =
    internal::Overlapping_region<Traits, Pair_range_3, First_of_pair_map, Second_of_pair_map, LSPF_region>;
    using LSPF_sorting =
    internal::Least_squares_plane_fit_sorting<Traits, Indices, Sphere_neighbor_query, Point_map_3>;
    using Region_growing_3_lspf =
    internal::Region_growing<Indices, Sphere_neighbor_query, LSPF_region, typename LSPF_sorting::Seed_map>;
    using Region_growing_3_olpf =
    internal::Region_growing<Indices, Sphere_neighbor_query, OLPF_region, typename LSPF_sorting::Seed_map>;

    using Groups = std::vector< std::set<std::size_t> >;

    Building_roofs_creator(
      const Point_map_3 point_map_3) :
    m_point_map_3(point_map_3)
    { }

    void create_cluster(
      const Indices& input,
      const FT region_growing_scale_3,
      const FT region_growing_angle_3,
      Indices& cluster) const {

      // Compute normals.
      std::vector<Vector_3> normals;
      Sphere_neighbor_query neighbor_query(
        input, region_growing_scale_3, m_point_map_3);
      Normal_estimator_3 estimator(
        input, neighbor_query, m_point_map_3);
      estimator.get_normals(normals);
      CGAL_assertion(normals.size() == input.size());

      // Remove vertical points.
      cluster.clear();
      const Vector_3 ref = Vector_3(FT(0), FT(0), FT(1));
      for (std::size_t i = 0; i < input.size(); ++i) {

        const auto& vec = normals[i];
        FT angle = angle_3d(vec, ref);
        if (angle > FT(90)) angle = FT(180) - angle;
        angle = FT(90) - angle;
        if (angle > region_growing_angle_3)
          cluster.push_back(input[i]);
      }

      std::vector<Point_3> roof_points;
      roof_points.reserve(cluster.size());
      for (const std::size_t idx : cluster)
        roof_points.push_back(get(m_point_map_3, idx));

      Saver<Traits> saver;
      saver.export_points(
        roof_points,
        Color(0, 0, 0),
        "/Users/monet/Documents/gf/lod/logs/buildings/tmp/roof-points");
    }

    void create_roof_regions(
      const Indices& cluster,
      const FT region_growing_scale_3,
      const FT region_growing_noise_level_3,
      const FT region_growing_angle_3,
      const FT region_growing_min_area_3,
      const FT region_growing_distance_to_line_3,
      const FT alpha_shape_size_2,
      const FT max_height_difference,
      std::vector<Indices>& roofs,
      Indices& unclassified) const {

      Sphere_neighbor_query neighbor_query(
        cluster, region_growing_scale_3, m_point_map_3);

      LSPF_sorting sorting(
        cluster, neighbor_query, m_point_map_3);
      sorting.sort();

      Pair_range_3 range;
      apply_default_region_growing(
        cluster, neighbor_query, sorting,
        region_growing_noise_level_3,
        region_growing_angle_3,
        region_growing_min_area_3,
        region_growing_distance_to_line_3,
        alpha_shape_size_2,
        roofs, unclassified, range);

      remove_overlapping_regions(
        cluster, neighbor_query, sorting, range,
        region_growing_noise_level_3,
        region_growing_angle_3,
        region_growing_min_area_3,
        region_growing_distance_to_line_3,
        alpha_shape_size_2,
        roofs, unclassified);

      const FT radius = region_growing_scale_3 * FT(2);
      const FT eps = max_height_difference / FT(2);
      std::map<std::size_t, bool> skip;
      reassign_indices(
        cluster, skip, radius, eps,
        unclassified, roofs);

      std::vector<Point_3> roof_points;
      std::vector< std::vector<Point_3> > roof_regions;
      roof_regions.reserve(roofs.size());
      for (const auto& roof : roofs) {
        roof_points.clear();
        for (const std::size_t idx : roof)
          roof_points.push_back(get(m_point_map_3, cluster[idx]));
        roof_regions.push_back(roof_points);
      }

      Saver<Traits> saver;
      saver.export_points(
        roof_regions,
        "/Users/monet/Documents/gf/lod/logs/buildings/tmp/roof-regions");
    }

  private:
    const Point_map_3 m_point_map_3;

    void reassign_indices(
      const Indices& cluster,
      const std::map<std::size_t, bool>& skip,
      const FT radius,
      const FT eps,
      Indices& indices,
      std::vector<Indices>& roofs) const {

      Sphere_neighbor_query neighbor_query(
        cluster, radius, m_point_map_3);

      Plane_3 plane;
      std::map<std::size_t, Plane_3> planes;
      for (std::size_t i = 0; i < roofs.size(); ++i) {
        if (skip.find(i) != skip.end()) continue;
        internal::plane_from_points_3(cluster, m_point_map_3, roofs[i], plane);
        planes[i] = plane;
      }

      std::map<std::size_t, std::size_t> plane_map;
      for (std::size_t i = 0; i < roofs.size(); ++i) {
        if (skip.find(i) != skip.end()) continue;
        for (const std::size_t idx : roofs[i])
          plane_map[idx] = i;
      }

      Indices neighbors;
      std::set<std::size_t> closest;
      std::vector<bool> states(indices.size(), false);

      for (std::size_t i = 0; i < indices.size(); ++i) {
        const std::size_t idx = indices[i];

        neighbors.clear();
        neighbor_query(idx, neighbors);
        closest.clear();
        for (const std::size_t neighbor : neighbors) {
          if (plane_map.find(neighbor) != plane_map.end())
            closest.insert(plane_map.at(neighbor));
        }
        if (closest.size() == 0) continue;

        FT min_dist = internal::max_value<FT>();
        std::size_t best_plane = std::size_t(-1);

        bool found = false;
        for (const std::size_t roof_idx : closest) {
          if (skip.find(roof_idx) != skip.end()) continue;

          const auto& plane = planes.at(roof_idx);
          const auto& p = get(m_point_map_3, cluster[idx]);
          const auto  q = plane.projection(p);

          const FT dist = internal::distance(p, q);
          if (dist < eps) {
            if (dist < min_dist) {
              min_dist = dist; best_plane = roof_idx;
            }
          } else { found = true; break; }
        }

        if (!found && best_plane != std::size_t(-1)) {
          roofs[best_plane].push_back(idx);
          states[i] = true;
        }
      }

      Indices updated;
      for (std::size_t i = 0; i < indices.size(); ++i)
        if (!states[i])
          updated.push_back(indices[i]);
      indices = updated;
    }

    void apply_default_region_growing(
      const Indices& cluster,
      Sphere_neighbor_query& neighbor_query,
      LSPF_sorting& sorting,
      const FT region_growing_noise_level_3,
      const FT region_growing_angle_3,
      const FT region_growing_min_area_3,
      const FT region_growing_distance_to_line_3,
      const FT alpha_shape_size_2,
      std::vector<Indices>& roofs,
      Indices& unclassified,
      Pair_range_3& range) const {

      std::vector<Vector_3> normals;
      Normal_estimator_3 estimator(
        cluster, neighbor_query, m_point_map_3);
      estimator.get_normals(normals);

      CGAL_assertion(cluster.size() == normals.size());

      range.clear();
      range.reserve(cluster.size());
      for (std::size_t i = 0; i < cluster.size(); ++i) {
        const Point_3& p = get(m_point_map_3, cluster[i]);
        range.push_back(std::make_pair(p, normals[i]));
      }

      First_of_pair_map point_map;
      Second_of_pair_map normal_map;
      LSPF_region lspf_region(
        range,
        region_growing_noise_level_3,
        region_growing_angle_3,
        region_growing_min_area_3,
        region_growing_distance_to_line_3,
        alpha_shape_size_2,
        point_map, normal_map);

      roofs.clear();
      Region_growing_3_lspf region_growing(
        cluster,
        neighbor_query,
        lspf_region,
        sorting.seed_map());
      region_growing.detect(std::back_inserter(roofs));

      unclassified.clear();
      region_growing.unassigned_items(std::back_inserter(unclassified));

      std::cout << "roofs before: " << roofs.size() << std::endl;
      std::cout << "unclassified before: " << unclassified.size() << std::endl;
    }

    void remove_overlapping_regions(
      const Indices& cluster,
      Sphere_neighbor_query& neighbor_query,
      LSPF_sorting& sorting,
      const Pair_range_3& range,
      const FT region_growing_noise_level_3,
      const FT region_growing_angle_3,
      const FT region_growing_min_area_3,
      const FT region_growing_distance_to_line_3,
      const FT alpha_shape_size_2,
      std::vector<Indices>& roofs,
      Indices& unclassified) const {

      First_of_pair_map point_map;
      Second_of_pair_map normal_map;
      LSPF_region lspf_region(
        range,
        region_growing_noise_level_3,
        region_growing_angle_3,
        region_growing_min_area_3,
        region_growing_distance_to_line_3,
        alpha_shape_size_2,
        point_map, normal_map);

      OLPF_region olpf_region(
        range, roofs, lspf_region,
        region_growing_noise_level_3,
        region_growing_angle_3,
        point_map, normal_map);

      std::vector<Indices> regions;
      Region_growing_3_olpf region_growing(
        cluster,
        neighbor_query,
        olpf_region,
        sorting.seed_map());

      region_growing.use_overlap(true);
      region_growing.detect(std::back_inserter(regions));

      Groups groups;
      olpf_region.get_groups(groups);
      fix_overlapping_points(range, groups, roofs, unclassified);

      /*
      std::map<std::size_t, bool> overlapping;
      olpf_region.get_overlapping_points(overlapping);
      remove_overlapping_points(overlapping, roofs, unclassified); */

      std::cout << "roofs after: " << roofs.size() << std::endl;
      std::cout << "unclassified after: " << unclassified.size() << std::endl;
    }

    void remove_overlapping_points(
      const std::map<std::size_t, bool>& overlapping,
      std::vector<Indices>& roofs,
      Indices& unclassified) const {

      std::vector<Indices> regions; Indices region;
      for (std::size_t k = 0; k < roofs.size(); ++k) {
        region.clear();
        for (const std::size_t idx : roofs[k]) {
          if (!overlapping.at(idx)) region.push_back(idx);
          else unclassified.push_back(idx);
        }
        if (region.size() > 0)
          regions.push_back(region);
      }
      roofs = regions;
    }

    void fix_overlapping_points(
      const Pair_range_3& range,
      const Groups& groups,
      std::vector<Indices>& roofs,
      Indices& unclassified) const {

      CGAL_assertion(range.size() == groups.size());
      std::vector<Indices> regions;
      create_regions(groups, roofs, regions);

      std::vector<Plane_3> planes;
      std::vector<Vector_3> normals;
      create_planes_and_normals(range, regions, planes, normals);

      relabel_points(range, groups, planes, normals, roofs, unclassified);
    }

    void create_regions(
      const Groups& groups,
      const std::vector<Indices>& roofs,
      std::vector<Indices>& regions) const {

      regions.clear();
      regions.resize(roofs.size());

      for (std::size_t i = 0; i < groups.size(); ++i) {
        const auto& group = groups[i];

        if (group.size() == 0 || group.size() > 1)
          continue;

        if (group.size() == 1) {
          const std::size_t k = *(group.begin());
          regions[k].push_back(i);
          continue;
        }
      }

      for (std::size_t k = 0; k < regions.size(); ++k) {
        auto& region = regions[k];
        if (region.size() < 3)
          region = roofs[k];
      }
    }

    void create_planes_and_normals(
      const Pair_range_3& range,
      const std::vector<Indices>& regions,
      std::vector<Plane_3>& planes,
      std::vector<Vector_3>& normals) const {

      planes.clear();
      planes.resize(regions.size());

      normals.clear();
      normals.resize(regions.size());

      for (std::size_t k = 0; k < regions.size(); ++k) {
        const auto& region = regions[k];

        Plane_3 plane; First_of_pair_map point_map_3;
        internal::plane_from_points_3(range, point_map_3, region, plane);
        planes[k] = plane;

        Vector_3 normal = plane.orthogonal_vector();
        const FT normal_length = static_cast<FT>(
        CGAL::sqrt(
          CGAL::to_double(
            normal.squared_length())));
        normal /= normal_length;
        normals[k] = normal;
      }
    }

    void relabel_points(
      const Pair_range_3& range,
      const Groups& groups,
      const std::vector<Plane_3>& planes,
      const std::vector<Vector_3>& normals,
      std::vector<Indices>& roofs,
      Indices& unclassified) const {

      std::vector<Indices> regions(roofs.size());
      for (std::size_t i = 0; i < groups.size(); ++i) {
        const auto& group = groups[i];

        if (group.size() == 0)
          continue;

        if (group.size() == 1) {
          const std::size_t k = *(group.begin());
          regions[k].push_back(i);
          continue;
        }

        if (group.size() > 1) {
          const auto& point = range[i].first;
          const auto& normal = range[i].second;
          const std::size_t k =
            estimate_best_roof_index(point, normal, group, planes, normals);
          if (k != std::size_t(-1))
            regions[k].push_back(i);
          else
            unclassified.push_back(i);
          continue;
        }
      }
      roofs = regions;
    }

    std::size_t estimate_best_roof_index(
      const Point_3& query_point,
      const Vector_3& query_normal,
      const std::set<std::size_t>& group,
      const std::vector<Plane_3>& planes,
      const std::vector<Vector_3>& normals) const {

      FT min_distance  =  internal::max_value<FT>();
      FT max_cos_value = -internal::max_value<FT>();

      std::size_t best_k = std::size_t(-1);
      for (const std::size_t k : group) {

        const auto& best_plane = planes[k];
        const auto& best_normal = normals[k];

        const FT distance = static_cast<FT>(
          CGAL::sqrt(
            CGAL::to_double(
              CGAL::squared_distance(query_point, best_plane))));

        const FT cos_value = CGAL::abs(
          CGAL::scalar_product(query_normal, best_normal));

        if (distance < min_distance && cos_value > max_cos_value) {
          min_distance = distance; max_cos_value = cos_value;
          best_k = k;
        }
      }
      return best_k;
    }
  };

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_INTERNAL_ROOFS_CREATOR_H
