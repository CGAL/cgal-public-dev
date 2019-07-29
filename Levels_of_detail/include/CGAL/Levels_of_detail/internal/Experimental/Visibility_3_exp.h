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

#ifndef CGAL_LEVELS_OF_DETAIL_INTERNAL_VISIBILITY_3_EXP_H
#define CGAL_LEVELS_OF_DETAIL_INTERNAL_VISIBILITY_3_EXP_H

// STL includes.
#include <vector>
#include <utility>

// CGAL includes.
#include <CGAL/assertions.h>

// Internal includes.
#include <CGAL/Levels_of_detail/enum.h>
#include <CGAL/Levels_of_detail/internal/struct.h>
#include <CGAL/Levels_of_detail/internal/utils.h>

// Spacial search.
#include <CGAL/Levels_of_detail/internal/Spacial_search/Nearest_face_neighbor_query.h>

// Shape detection.
#include <CGAL/Levels_of_detail/internal/Shape_detection/Region_growing.h>
#include <CGAL/Levels_of_detail/internal/Shape_detection/Coplanar_region.h>

// Testing.
#include "../../../../../test/Levels_of_detail/include/Saver.h"

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

  template<
  typename GeomTraits,
  typename InputRange,
  typename PointMap>
  class Visibility_3_exp {
			
  public:
    using Traits = GeomTraits;
    using Input_range = InputRange;
    using Point_map = PointMap;

    using FT = typename Traits::FT;
    using Point_2 = typename Traits::Point_2;
    using Point_3 = typename Traits::Point_3;
    using Vector_3 = typename Traits::Vector_3;
    using Line_3 = typename Traits::Line_3;
    using Plane_3 = typename Traits::Plane_3;
    using Triangle_2 = typename Traits::Triangle_2;

    using Indices = std::vector<std::size_t>;
    using Partition_3 = internal::Partition_3<Traits>;
    using Edge = typename Partition_3::Edge;
    using Building = internal::Building<Traits>;
    using Points_3 = std::vector<Point_3>;
    using Polygons = std::vector<Points_3>;

    using Saver = Saver<Traits>;
    using Color = CGAL::Color;

    using Nearest_face_neighbor_query = internal::Nearest_face_neighbor_query<Traits>;
    using Coplanar_region = internal::Coplanar_region<Traits>;
    using Region_growing = internal::Region_growing<Polygons, Nearest_face_neighbor_query, Coplanar_region>;

    Visibility_3_exp(
      const Input_range& input_range,
      const Point_map& point_map,
      const Building& building,
      const std::vector<Indices>& roof_points_3) :
    m_input_range(input_range),
    m_point_map(point_map),
    m_building(building),
    m_roof_points_3(roof_points_3)
    { }

    void compute(Partition_3& partition) const {
      
      std::vector<Points_3> polygons;
      const std::vector<Edge>& edges = partition.edges;
      remove_unnecessary_edges(edges, polygons);
      
      std::vector<Polygons> groups;
      group_polygons(polygons, groups);

      std::vector<int> roofs(groups.size(), -1);
      assign_roof_indices(groups, roofs);

      Saver saver;
      std::vector<Point_3> points;
      for (std::size_t i = 0; i < groups.size(); ++i) {

        Random rand(i);
        const auto r = static_cast<unsigned char>(64 + rand.get_int(0, 192));
        const auto g = static_cast<unsigned char>(64 + rand.get_int(0, 192));
        const auto b = static_cast<unsigned char>(64 + rand.get_int(0, 192));
        const Color color = Color(r, g, b);

        points.clear();
        const int roof_idx = roofs[i];
        if (roof_idx >= 0) {
          for (const std::size_t idx : m_roof_points_3[roof_idx]) {
            const auto& p = get(m_point_map, *(m_input_range.begin() + idx));
            points.push_back(p);
          }
        }

        const std::string rname = 
        "/Users/monet/Documents/lod/logs/buildings/groups/polygon-soup-" + std::to_string(i);
        saver.export_polygon_soup(groups[i], color, rname);
        const std::string pname =
        "/Users/monet/Documents/lod/logs/buildings/groups/point-set-" + std::to_string(i);
        saver.export_points(points, color, pname);
      }
    }
    
  private:
    const Input_range& m_input_range;
    const Point_map& m_point_map;
    const Building& m_building;
    const std::vector<Indices>& m_roof_points_3;
    CGAL::Bbox_3 m_bbox;

    void remove_unnecessary_edges(
      const std::vector<Edge>& edges,
      std::vector<Points_3>& polygons) const {
      
      for (const auto& edge : edges) {
        if (is_vertical(edge)) continue;
        if (is_below_building(edge)) continue;
        if (is_above_building(edge)) continue;
        if (is_outside_boundary(edge)) continue;

        polygons.push_back(edge.polygon);
      }
    }

    bool is_vertical(const Edge& edge) const {
      Vector_3 norm;
      const Vector_3 vert = Vector_3(FT(0), FT(0), FT(1));
      internal::compute_normal_3(edge.polygon, norm);
      const FT angle = internal::angle_3d(vert, norm);
      return (angle >= FT(89) && angle <= FT(91));
    }

    bool is_below_building(const Edge& edge) const {
      Point_3 b;
      internal::compute_barycenter_3(edge.polygon, b);
      return b.z() < m_building.bottom_z;
    }

    bool is_above_building(const Edge& edge) const {
      Point_3 b;
      internal::compute_barycenter_3(edge.polygon, b);
      return b.z() > m_building.top_z;
    }

    bool is_outside_boundary(const Edge& edge) const {

      Point_3 b;
      internal::compute_barycenter_3(edge.polygon, b);
      const Point_2 p = internal::point_2_from_point_3(b);
      const auto& tri = m_building.base1.triangulation.delaunay;
      const auto fh = tri.locate(p);
      return !fh->info().tagged;
    }

    void group_polygons(
      const Polygons& polygons,
      std::vector<Polygons>& groups) const {

      Nearest_face_neighbor_query neighbor_query(polygons);
      Coplanar_region region(polygons);
      Region_growing region_growing(
        polygons, neighbor_query, region);
      std::vector<Indices> regions;
      region_growing.detect(std::back_inserter(regions));

      groups.resize(regions.size());
      for (std::size_t i = 0; i < regions.size(); ++i) {
        for (const std::size_t idx : regions[i])
          groups[i].push_back(polygons[idx]);
      }
    }

    void assign_roof_indices(
      const std::vector<Polygons>& groups,
      std::vector<int>& roofs) const {

      std::vector<Vector_3> norms;

      Plane_3 plane;
      for (const auto& region : m_roof_points_3) {
        internal::plane_from_points_3(
          m_input_range, m_point_map, region, plane);
        
        Vector_3 norm = plane.orthogonal_vector();
        internal::normalize(norm);
        norms.push_back(norm);
      }

      
    }
  };

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_INTERNAL_VISIBILITY_3_EXP_H
