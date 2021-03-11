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

#ifndef CGAL_LEVELS_OF_DETAIL_INTERNAL_VISIBILITY_3_EXP_1_H
#define CGAL_LEVELS_OF_DETAIL_INTERNAL_VISIBILITY_3_EXP_1_H

// Example:

// void compute_visibility_3_exp_1(const FT) {

//   using Visibility_3_exp_1 = internal::Visibility_3_exp_1<
//   Traits, Points_3, Point_map_3>;

//   if (m_partition_3.empty()) return;
//   Visibility_3_exp_1 visibility(
//     m_cluster,
//     m_data.point_map_3,
//     m_building,
//     m_roof_points_3,
//     m_data.parameters.buildings.alpha_shape_size_2,
//     m_data.parameters.buildings.graphcut_beta_3);
//   visibility.compute(m_partition_3);

//   std::cout << "visibility finished" << std::endl;
// }

// STL includes.
#include <vector>
#include <utility>

// CGAL includes.
#include <CGAL/assertions.h>

// Internal includes.
#include <CGAL/Levels_of_detail/enum.h>
#include <CGAL/Levels_of_detail/internal/struct.h>
#include <CGAL/Levels_of_detail/internal/utils.h>

// Visibility.
#include <CGAL/Levels_of_detail/internal/Visibility/Visibility_2.h>

// Graphcut.
#include <CGAL/Levels_of_detail/internal/Graphcut/Graphcut.h>

// Testing.
#include "../../../../../test/Levels_of_detail/include/Saver.h"

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

  template<
  typename GeomTraits,
  typename InputRange,
  typename PointMap>
  class Visibility_3_exp_1 {

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
    using Triangle_3 = typename Traits::Triangle_3;

    using Indices = std::vector<std::size_t>;
    using Partition_2 = internal::Partition_2<Traits>;
    using Partition_3 = internal::Partition_3<Traits>;

    using Edge = typename Partition_3::Edge;
    using Face = typename Partition_2::Face;

    using Building = internal::Building<Traits>;
    using Roof = typename Building::Roof;
    using Wall = typename Building::Wall;

    using Points_2 = std::vector<Point_2>;
    using Points_3 = std::vector<Point_3>;

    using Identity_map_2 = CGAL::Identity_property_map<Point_2>;
    using Identity_map_3 = CGAL::Identity_property_map<Point_3>;

    using Point_map_2 =
    internal::Item_property_map<Points_2, Identity_map_2>;
    using Visibility_2 = internal::Visibility_2<Traits, Point_map_2>;
    using Graphcut_2 = internal::Graphcut<Traits, Partition_2>;

    using Saver = Saver<Traits>;
    using Color = CGAL::Color;

    Visibility_3_exp_1(
      const Input_range& input_range,
      const Point_map& point_map,
      Building& building, // temporary solution
      const std::vector<Indices>& roof_points_3,
      const FT alpha_shape_size_2,
      const FT graphcut_beta_2) :
    m_input_range(input_range),
    m_point_map(point_map),
    m_building(building),
    m_roof_points_3(roof_points_3),
    m_alpha_shape_size_2(alpha_shape_size_2),
    m_graphcut_beta_2(graphcut_beta_2)
    { }

    void compute(Partition_3& partition) {

      m_building.walls2.clear();
      m_building.roofs2.clear();

      std::vector<Indices> groups;
      std::vector<Points_3> polygons;

      const std::vector<Edge>& edges = partition.edges;
      remove_unnecessary_edges(edges, groups, polygons);
      apply_2d_visibility(polygons, groups);

      Points_3 points;
      std::vector<Points_3> group;
      for (std::size_t i = 0; i < groups.size(); ++i) {

        Random rand(i);
        const auto r = static_cast<unsigned char>(64 + rand.get_int(0, 192));
        const auto g = static_cast<unsigned char>(64 + rand.get_int(0, 192));
        const auto b = static_cast<unsigned char>(64 + rand.get_int(0, 192));
        const Color color = Color(r, g, b);

        group.clear();
        for (const std::size_t idx : groups[i])
          group.push_back(polygons[idx]);

        points.clear();
        for (const std::size_t idx : m_roof_points_3[i]) {
          const auto& p = get(m_point_map, *(m_input_range.begin() + idx));
          points.push_back(p);
        }

        // Saver saver;
        // const std::string rname =
        // "/Users/monet/Documents/lod/logs/buildings/groups/polygon-soup-" + std::to_string(i);
        // saver.export_polygon_soup(group, color, rname);
        // const std::string pname =
        // "/Users/monet/Documents/lod/logs/buildings/groups/point-set-" + std::to_string(i);
        // saver.export_points(points, color, pname);

        add_walls_and_roofs(group); // hack!
      }
    }

  private:
    const Input_range& m_input_range;
    const Point_map& m_point_map;
    Building& m_building;
    const std::vector<Indices>& m_roof_points_3;
    const FT m_alpha_shape_size_2;
    const FT m_graphcut_beta_2;

    void add_walls_and_roofs(const std::vector<Points_3>& polygons) {

      auto& roofs = m_building.roofs2;
      auto& walls = m_building.walls2;

      // Base.
      m_building.base2 = m_building.base1;

      Roof roof; Wall wall;
      for (const auto& polygon : polygons) {
        const Point_3& ref = polygon[0];

        // Roofs.
        roof.triangles.clear();
        for (std::size_t i = 1; i < polygon.size() - 1; ++i) {
          const std::size_t ip = i + 1;

          const Point_3& p1 = ref;
          const Point_3& p2 = polygon[i];
          const Point_3& p3 = polygon[ip];
          const Triangle_3 triangle = Triangle_3(p1, p2, p3);
          roof.triangles.push_back(triangle);
        }
        roofs.push_back(roof);

        // Walls.
        wall.triangles.clear();
        for (std::size_t i = 0; i < polygon.size(); ++i) {
          const std::size_t ip = (i + 1) % polygon.size();

          const auto& p1 = polygon[i];
          const auto& p2 = polygon[ip];

          const Triangle_3 tri1 = Triangle_3(
            Point_3(p1.x(), p1.y(), m_building.bottom_z),
            Point_3(p2.x(), p2.y(), m_building.bottom_z),
            p2);
          const Triangle_3 tri2 = Triangle_3(
            p2,
            p1,
            Point_3(p1.x(), p1.y(), m_building.bottom_z));
          wall.triangles.push_back(tri1);
          wall.triangles.push_back(tri2);
        }
        walls.push_back(wall);
      }
    }

    void remove_unnecessary_edges(
      const std::vector<Edge>& edges,
      std::vector<Indices>& groups,
      std::vector<Points_3>& polygons) const {

      groups.resize(edges.back().plane_index + 1);
      Point_3 b; std::size_t poly_idx = 0;
      for (const auto& edge : edges) {

        if (edge.plane_index < 0) continue;
        internal::compute_barycenter_3(edge.polygon, b);
        if (is_below_building(b, edge)) continue;
        if (is_above_building(b, edge)) continue;
        if (is_outside_boundary(b, edge)) continue;

        polygons.push_back(edge.polygon);
        groups[edge.plane_index].push_back(poly_idx);
        ++poly_idx;
      }

      std::vector<Indices> new_groups;
      for (const auto& group : groups)
        if (!group.empty())
          new_groups.push_back(group);
      groups = new_groups;
    }

    bool is_below_building(const Point_3& b, const Edge& edge) const {
      return b.z() < m_building.bottom_z;
    }

    bool is_above_building(const Point_3& b, const Edge& edge) const {
      return b.z() > m_building.top_z;
    }

    bool is_outside_boundary(const Point_3& b, const Edge& edge) const {
      const Point_2 p = internal::point_2_from_point_3(b);
      const auto& tri = m_building.base1.triangulation.delaunay;
      const auto fh = tri.locate(p);
      return !fh->info().tagged;
    }

    void apply_2d_visibility(
      const std::vector<Points_3>& polygons,
      std::vector<Indices>& groups) const {

      Points_3 tmp_points;
      std::vector<Points_3> tmp_polygons;
      for (std::size_t i = 0; i < groups.size(); ++i) {

        tmp_polygons.clear();
        for (const std::size_t idx : groups[i])
          tmp_polygons.push_back(polygons[idx]);

        tmp_points.clear();
        for (const std::size_t idx : m_roof_points_3[i]) {
          const auto& p = get(m_point_map, *(m_input_range.begin() + idx));
          tmp_points.push_back(p);
        }

        apply_region_2d_visibility(i, tmp_points, tmp_polygons, groups[i]);
      }
    }

    void apply_region_2d_visibility(
      const std::size_t id,
      Points_3& points,
      std::vector<Points_3>& polygons,
      Indices& group) const {

      Point_3 b;
      compute_barycenter(points, polygons, b);
      Plane_3 plane;
      Identity_map_3 identity_map_3;
      internal::plane_from_points_3(points, identity_map_3, plane);

      Points_2 points_2d; Points_3 points_save;
      create_2d_points(points, b, plane, points_2d, points_save);

      std::vector<Points_2> polygons_2d;
      std::vector<Points_3> polygons_save;
      create_2d_polygons(polygons, b, plane, polygons_2d, polygons_save);

      // Saver saver;
      // const Color color = Color(0, 0, 0);
      // const std::string rname =
      // "/Users/monet/Documents/lod/logs/buildings/rotated/polygon-soup-" + std::to_string(id);
      // saver.export_polygon_soup(polygons_save, color, rname);
      // const std::string pname =
      // "/Users/monet/Documents/lod/logs/buildings/rotated/point-set-" + std::to_string(id);
      // saver.export_points(points_save, color, pname);

      Partition_2 partition_2;
      compute_visibility_2(id, points_2d, polygons_2d, partition_2);

      Indices new_group;
      for (std::size_t i = 0; i < partition_2.faces.size(); ++i)
        if (partition_2.faces[i].visibility == Visibility_label::INSIDE)
          new_group.push_back(group[i]);
      group = new_group;
    }

    void compute_barycenter(
      const Points_3& points,
      const std::vector<Points_3>& polygons,
      Point_3& b) const {

      Point_3 b1;
      internal::compute_barycenter_3(points, b1);

      Point_3 b2;
      Indices indices;
      indices.reserve(polygons.size());
      for (std::size_t i = 0; i < polygons.size(); ++i)
        indices.push_back(i);
      internal::compute_barycenter_3(polygons, indices, b2);

      const FT x = (b1.x() + b2.x()) / FT(2);
      const FT y = (b1.y() + b2.y()) / FT(2);
      const FT z = (b1.z() + b2.z()) / FT(2);

      b = Point_3(x, y, z);
    }

    void create_2d_points(
      const Points_3& points,
      const Point_3& b,
      const Plane_3& plane,
      Points_2& points_2d,
      Points_3& points_save) const {

      for (const auto& p : points) {
        const Point_2 q = internal::to_2d(p, b, plane);
        points_2d.push_back(q);
        points_save.push_back(Point_3(q.x(), q.y(), FT(0)));
      }
    }

    void create_2d_polygons(
      const std::vector<Points_3>& polygons,
      const Point_3& b,
      const Plane_3& plane,
      std::vector<Points_2>& polygons_2d,
      std::vector<Points_3>& polygons_save) const {

      polygons_2d.resize(polygons.size());
      polygons_save.resize(polygons.size());

      for (std::size_t i = 0; i < polygons.size(); ++i) {
        for (const auto& p : polygons[i]) {
          const Point_2 q = internal::to_2d(p, b, plane);
          polygons_2d[i].push_back(q);
          polygons_save[i].push_back(Point_3(q.x(), q.y(), FT(0)));
        }
      }
    }

    void compute_visibility_2(
      const std::size_t id,
      const Points_2& points,
      const std::vector<Points_2>& polygons,
      Partition_2& partition_2) const {

      Indices stub, indices;
      for (std::size_t i = 0; i < points.size(); ++i)
        indices.push_back(i);
      Identity_map_2 identity_map_2;
      Point_map_2 point_map_2(points, identity_map_2);

      const Visibility_2 visibility(
        stub, indices,
        point_map_2,
        m_alpha_shape_size_2, FT(1) / FT(5));

      for (const auto& polygon : polygons)
        partition_2.faces.push_back(Face(polygon));
      visibility.compute(partition_2);

      const Graphcut_2 graphcut(m_graphcut_beta_2);
      graphcut.apply(partition_2);
    }
  };

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_INTERNAL_VISIBILITY_3_EXP_1_H
