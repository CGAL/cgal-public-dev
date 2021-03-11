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

#ifndef CGAL_LEVELS_OF_DETAIL_INTERNAL_WALLS_CREATOR_H
#define CGAL_LEVELS_OF_DETAIL_INTERNAL_WALLS_CREATOR_H

// STL includes.
#include <vector>

// CGAL includes.
#include <CGAL/assertions.h>
#include <CGAL/property_map.h>
#include <CGAL/Optimal_transportation_reconstruction_2.h>

// Internal includes.
#include <CGAL/Levels_of_detail/internal/utils.h>
#include <CGAL/Levels_of_detail/internal/struct.h>

// Spatial search.
#include <CGAL/Levels_of_detail/internal/Spatial_search/K_neighbor_query.h>
#include <CGAL/Levels_of_detail/internal/Spatial_search/Sphere_neighbor_query.h>

// Shape detection.
#include <CGAL/Levels_of_detail/internal/Shape_detection/Region_growing.h>
#include <CGAL/Levels_of_detail/internal/Shape_detection/Estimate_normals_2.h>
#include <CGAL/Levels_of_detail/internal/Shape_detection/Least_squares_line_fit_region.h>
#include <CGAL/Levels_of_detail/internal/Shape_detection/Least_squares_line_fit_sorting.h>

// Testing.
#include "../../../../../test/Levels_of_detail/include/Saver.h"

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

  template<
  typename GeomTraits>
  class Building_walls_creator {

  public:
    using Traits = GeomTraits;

    using FT = typename Traits::FT;
    using Point_2 = typename Traits::Point_2;
    using Point_3 = typename Traits::Point_3;
    using Vector_2 = typename Traits::Vector_2;
    using Vector_3 = typename Traits::Vector_3;
    using Segment_2 = typename Traits::Segment_2;
    using Line_2 = typename Traits::Line_2;

    using Points_2 = std::vector<Point_2>;
    using Vectors_2 = std::vector<Vector_2>;

    using Indices = std::vector<std::size_t>;
    using Identity_map = CGAL::Identity_property_map<Point_2>;

    using Points_3 = std::vector<Point_3>;
    using Vectors_3 = std::vector<Vector_3>;

    using Pair_item_2 = std::pair<Point_2, Vector_2>;
    using Pair_range_2 = std::vector<Pair_item_2>;
    using First_of_pair_map = CGAL::First_of_pair_property_map<Pair_item_2>;
    using Second_of_pair_map = CGAL::Second_of_pair_property_map<Pair_item_2>;

    using Boundary_point_map_2 =
    internal::Item_property_map<Points_2, Identity_map>;

    using K_neighbor_query =
    internal::K_neighbor_query<Traits, Points_2, Identity_map>;
    using Sphere_neighbor_query =
    internal::Sphere_neighbor_query<Traits, Points_2, Identity_map>;

    using Saver = Saver<Traits>;
    using Color = CGAL::Color;

    using Otr = CGAL::Optimal_transportation_reconstruction_2<Traits, Identity_map>;

    Building_walls_creator(
      const Points_2& boundary_points_2) :
    m_boundary_points_2(boundary_points_2)
    { }

    void create_wall_and_roof_points(
      const FT region_growing_scale_3,
      const FT region_growing_angle_3,
      const Points_3& input_points,
      Points_3& wall_points, Points_3& roof_points,
      Vectors_3& wall_normals, Vectors_3& roof_normals) {

      using Identity_map_3 =
        CGAL::Identity_property_map<Point_3>;
      using SNQ =
        internal::Sphere_neighbor_query<Traits, Points_3, Identity_map_3>;
      using NE3 =
        internal::Estimate_normals_3<Traits, Points_3, Identity_map_3, SNQ>;

      // Compute normals.
      Identity_map_3 identity_map_3;
      Vectors_3 input_normals;
      SNQ neighbor_query(
        input_points, region_growing_scale_3, identity_map_3);
      NE3 estimator(
        input_points, neighbor_query, identity_map_3);
      estimator.get_normals(input_normals);
      CGAL_assertion(input_normals.size() == input_points.size());

      // Create wall and roof points.
      create_wall_and_roof_points(
        region_growing_angle_3,
        input_points, input_normals,
        wall_points, roof_points,
        wall_normals, roof_normals);
    }

    void create_wall_and_roof_points(
      const FT region_growing_angle_3,
      const Points_3& input_points,
      const Vectors_3& input_normals,
      Points_3& wall_points, Points_3& roof_points,
      Vectors_3& wall_normals, Vectors_3& roof_normals) {

      wall_points.clear(); roof_points.clear();
      wall_normals.clear(); roof_normals.clear();

      const Vector_3 ref = Vector_3(FT(0), FT(0), FT(1));
      for (std::size_t i = 0; i < input_points.size(); ++i) {

        const auto& p = input_points[i];
        const auto& n = input_normals[i];

        FT angle = angle_3d(n, ref);
        if (angle > FT(90)) angle = FT(180) - angle;
        angle = FT(90) - angle;
        if (angle <= region_growing_angle_3) {

          wall_points.push_back(p);
          wall_normals.push_back(n);

        } else {

          roof_points.push_back(p);
          roof_normals.push_back(n);
        }
      }
    }

    void create_wall_regions(
      const FT region_growing_scale_2,
      const FT region_growing_noise_level_2,
      const FT region_growing_angle_2,
      const FT region_growing_min_length_2,
      std::vector<Indices>& regions) {

      Identity_map identity_map;
      Sphere_neighbor_query neighbor_query(
        m_boundary_points_2, region_growing_scale_2, identity_map);
      apply_region_growing(
        region_growing_noise_level_2,
        region_growing_angle_2,
        region_growing_min_length_2,
        neighbor_query,
        regions);
    }

    template<typename Neighbor_query>
    void create_wall_regions(
      const FT region_growing_noise_level_2,
      const FT region_growing_angle_2,
      const FT region_growing_min_length_2,
      Neighbor_query& neighbor_query,
      std::vector<Indices>& regions) {

      apply_region_growing(
        region_growing_noise_level_2,
        region_growing_angle_2,
        region_growing_min_length_2,
        neighbor_query,
        regions);
    }

    void create_boundaries(
      const std::vector<Indices>& regions,
      std::vector<Segment_2>& segments) {

      segments.clear();
      segments.reserve(regions.size());

      Identity_map identity_map;
      Boundary_point_map_2 point_map_2(
        m_boundary_points_2, identity_map);

      Line_2 line; Point_2 p, q;
      std::vector<std::size_t> indices;
      for (const auto& item_range : regions) {
        indices.clear();
        for (std::size_t i = 0; i < item_range.size(); ++i)
          indices.push_back(i);
        internal::line_from_points_2(
          item_range, point_map_2, line);
        internal::boundary_points_on_line_2(
          item_range, point_map_2, indices, line, p, q);
        segments.push_back(Segment_2(p, q));
      }
      CGAL_assertion(
        segments.size() == regions.size());
    }

    void create_boundaries(
      const std::vector<Indices>& regions,
      std::vector<Segment_2>& segments,
      std::vector< std::vector<Point_2> >& projected) {

      segments.clear();
      segments.reserve(regions.size());

      projected.clear();
      projected.resize(regions.size());

      Identity_map identity_map;
      Boundary_point_map_2 point_map_2(
        m_boundary_points_2, identity_map);

      Line_2 line; Point_2 p, q;
      std::vector<std::size_t> indices;
      for (std::size_t j = 0; j < regions.size(); ++j) {
        const auto& item_range = regions[j];

        indices.clear();
        for (std::size_t i = 0; i < item_range.size(); ++i)
          indices.push_back(i);
        internal::line_from_points_2(
          item_range, point_map_2, line);
        internal::boundary_points_on_line_2(
          item_range, point_map_2, indices, line, p, q);
        internal::project_on_line_2(
          item_range, point_map_2, indices, line, projected[j]);
        segments.push_back(Segment_2(p, q));
      }
      CGAL_assertion(
        segments.size() == regions.size());
    }

    void save_polylines(
      const std::vector<Segment_2>& segments,
      const std::string name) {

      CGAL_assertion(segments.size() > 0);
      std::vector< std::vector<Point_3> > polylines(segments.size());
      for (std::size_t i = 0; i < segments.size(); ++i) {
        const Point_2& s = segments[i].source();
        const Point_2& t = segments[i].target();

        polylines[i].push_back(Point_3(s.x(), s.y(), FT(0)));
        polylines[i].push_back(Point_3(t.x(), t.y(), FT(0)));
      }

      Saver saver;
      saver.export_polylines(polylines, name);
    }

    void reconstruct_with_optimal_transport(
      const FT noise_level,
      std::vector<Segment_2>& segments) {

      if (m_boundary_points_2.empty())
        return;

      segments.clear();

      Identity_map identity_map;
      Otr otr(m_boundary_points_2, identity_map);
      otr.run_under_wasserstein_tolerance(noise_level);
      otr.list_output(
        boost::make_function_output_iterator([&](const Point_2&) -> void { }),
        boost::make_function_output_iterator([&](const Segment_2& segment) -> void {
          segments.push_back(segment);
        })
      );
    }

  private:
    const Points_2& m_boundary_points_2;

    template<typename Neighbor_query>
    void apply_region_growing(
      const FT region_growing_noise_level_2,
      const FT region_growing_angle_2,
      const FT region_growing_min_length_2,
      Neighbor_query& neighbor_query,
      std::vector<Indices>& regions) {

      regions.clear();
      Identity_map identity_map;

      using Normal_estimator_2 =
      internal::Estimate_normals_2<Traits, Points_2, Identity_map, Neighbor_query>;
      using LSLF_region =
      internal::Least_squares_line_fit_region<Traits, Pair_range_2, First_of_pair_map, Second_of_pair_map>;
      using LSLF_sorting =
      internal::Least_squares_line_fit_sorting<Traits, Points_2, Neighbor_query, Identity_map>;
      using Region_growing_2 =
      internal::Region_growing<Points_2, Neighbor_query, LSLF_region, typename LSLF_sorting::Seed_map>;

      Vectors_2 normals;
      Normal_estimator_2 estimator(
        m_boundary_points_2, neighbor_query, identity_map);
      estimator.get_normals(normals);

      CGAL_assertion(m_boundary_points_2.size() == normals.size());
      Pair_range_2 range;
      range.reserve(m_boundary_points_2.size());
      for (std::size_t i = 0; i < m_boundary_points_2.size(); ++i)
        range.push_back(std::make_pair(m_boundary_points_2[i], normals[i]));

      First_of_pair_map point_map;
      Second_of_pair_map normal_map;
      LSLF_region region(
        range,
        region_growing_noise_level_2,
        region_growing_angle_2,
        region_growing_min_length_2,
        point_map,
        normal_map);

      LSLF_sorting sorting(
        m_boundary_points_2, neighbor_query, identity_map);
      sorting.sort();

      Region_growing_2 region_growing(
        m_boundary_points_2,
        neighbor_query,
        region,
        sorting.seed_map());
      region_growing.detect(std::back_inserter(regions));
    }
  };

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_INTERNAL_WALLS_CREATOR_H
