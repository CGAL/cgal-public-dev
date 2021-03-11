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

#ifndef CGAL_LEVELS_OF_DETAIL_INTERNAL_VISIBILITY_3_DEPR_H
#define CGAL_LEVELS_OF_DETAIL_INTERNAL_VISIBILITY_3_DEPR_H

// STL includes.
#include <vector>
#include <utility>

// CGAL includes.
#include <CGAL/assertions.h>

// Internal includes.
#include <CGAL/Levels_of_detail/enum.h>
#include <CGAL/Levels_of_detail/internal/struct.h>
#include <CGAL/Levels_of_detail/internal/utils.h>

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

  template<
  typename GeomTraits,
  typename InputRange,
  typename PointMap>
  class Visibility_3_depr {

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
    using Stats = std::pair<FT, FT>;
    using Face = typename Partition_3::Face;
    using Building = internal::Building<Traits>;

    Visibility_3_depr(
      const Input_range& input_range,
      const Point_map& point_map,
      const Building& building,
      const std::vector<Indices>& roof_points_3,
      const FT visibility_scale_3) :
    m_input_range(input_range),
    m_point_map(point_map),
    m_building(building),
    m_roof_points_3(roof_points_3),
    m_distance_tolerance(visibility_scale_3),
    m_angle_threshold(FT(10)),
    m_height_offset(m_distance_tolerance / FT(20))
    { }

    void compute(Partition_3& partition) const {

      if (partition.empty()) return;
      label_exterior_faces(partition.faces);
      for (auto& face : partition.faces) {
        if (face.exterior) {
          face.visibility = Visibility_label::OUTSIDE;
          face.inside = FT(0);
          face.outside = FT(1);
        } else compute_face_label(face);
      }
    }

  private:
    const Input_range& m_input_range;
    const Point_map& m_point_map;
    const Building& m_building;
    const std::vector<Indices>& m_roof_points_3;
    const FT m_distance_tolerance;

    // Internal parameters.
    const FT m_angle_threshold;
    const FT m_height_offset;

    void label_exterior_faces(
      std::vector<Face>& faces) const {

      for (auto& face : faces) {
        face.exterior = false;
        const auto& neighbors = face.neighbors;
        for (const int idx : neighbors) {
          if (idx < 0) {
            face.exterior = true;
            break;
          }
        }
      }
    }

    void compute_face_label(Face& face) const {

      const Stats stats = estimate_in_out_values(face);
      CGAL_assertion(stats.first  >= FT(0) && stats.first  <= FT(1));
      CGAL_assertion(stats.second >= FT(0) && stats.second <= FT(1));
      CGAL_assertion(
        CGAL::abs(stats.first + stats.second - FT(1)) < internal::tolerance<FT>());

      if (stats.first > FT(1) / FT(2))
        face.visibility = Visibility_label::INSIDE;
      else
        face.visibility = Visibility_label::OUTSIDE;
      face.inside = stats.first;
      face.outside = stats.second;
    }

    Stats estimate_in_out_values(const Face& polyhedron) const {

      Point_3 b;
      internal::compute_barycenter_3(polyhedron.vertices, b);
      if (is_above_building(b))
        return std::make_pair(FT(0), FT(1));
      if (is_below_building(b))
        return std::make_pair(FT(0), FT(1));
      if (is_out_of_building(b))
        return std::make_pair(FT(1) / FT(5), FT(4) / FT(5));
      if (has_vertices_outside(polyhedron))
        return std::make_pair(FT(2) / FT(5), FT(3) / FT(5));
      return estimate_in_out_values_statistically(polyhedron, b);
    }

    bool is_above_building(const Point_3& query) const {
      return query.z() > m_building.top_z;
    }

    bool is_below_building(const Point_3& query) const {
      return query.z() < m_building.bottom_z;
    }

    bool is_out_of_building(const Point_3& query) const {

      const Point_2 p = Point_2(query.x(), query.y());

      const auto& tri = m_building.base1.triangulation.delaunay;
      for (auto fh = tri.finite_faces_begin();
      fh != tri.finite_faces_end(); ++fh) {
        if (!fh->info().interior) continue;
        const Triangle_2 triangle = Triangle_2(
          fh->vertex(0)->point(),
          fh->vertex(1)->point(),
          fh->vertex(2)->point());

        if (internal::is_within_triangle_2(p, triangle, m_distance_tolerance))
          return false;
      }
      return true;
    }

    bool has_vertices_outside(const Face& polyhedron) const {

      std::size_t count = 0;
      for (const auto& p : polyhedron.vertices) {
        const bool is_out = is_out_of_building(p);
        if (is_out) ++count;
        if (is_out && count > 0)
          return true;
      }
      return false;
    }

    Stats estimate_in_out_values_statistically(
      const Face& polyhedron, const Point_3& b) const {

      std::size_t in = 0, out = 0;
      std::vector<std::size_t> indices;

      std::vector<Point_3> polygon;
      for (const auto& face : polyhedron.faces) {

        polygon.clear();
        for (const std::size_t idx : face)
          polygon.push_back(polyhedron.vertices[idx]);
        if (internal::is_vertical_polygon(polygon, m_angle_threshold))
          continue;
        process_face(polygon, indices, in, out);
      }
      process_middle_plane(b, indices, in, out);

      if (in == 0 && out == 0)
        return std::make_pair(FT(1) / FT(5), FT(4) / FT(5));

      const FT tmp_in = static_cast<FT>(in);
      const FT tmp_out = static_cast<FT>(out);
      const FT sum = tmp_in + tmp_out;
      CGAL_assertion(sum > FT(0));

      const FT final_in = tmp_in  / sum;
      const FT final_out = tmp_out / sum;
      return std::make_pair(final_in, final_out);
    }

    void process_face(
      const std::vector<Point_3>& poly_3,
      std::vector<std::size_t>& indices,
      std::size_t& in,
      std::size_t& out) const {

      std::vector<Point_2> poly_2;
      internal::polygon_3_to_polygon_2(poly_3, poly_2);
      if (!CGAL::is_simple_2(poly_2.begin(), poly_2.end()))
        return;

      for (const auto& region : m_roof_points_3) {
        for (const std::size_t idx : region) {
          const auto& p = get(m_point_map, *(m_input_range.begin() + idx));

          FT z = internal::max_value<FT>();
          const Point_2 query = Point_2(p.x(), p.y());
          if (internal::is_inside_polygon_2(query, poly_2))
            z = internal::intersect_with_polygon_3(p, poly_3,
            m_building.bottom_z, m_building.top_z);
          if (z == internal::max_value<FT>())
            continue;

          indices.push_back(idx);
          if (is_inside_building(z, p.z())) ++in;
          else ++out;
        }
      }
    }

    void process_middle_plane(
      const Point_3& b,
      const std::vector<std::size_t>& indices,
      size_t& in,
      size_t& out) const {

      Line_3 line;
      Plane_3 middle_plane(b, Vector_3(FT(0), FT(0), FT(1)));
      for (const std::size_t idx : indices) {
        const auto& p = get(m_point_map, *(m_input_range.begin() + idx));

        const Point_3 p1 = Point_3(p.x(), p.y(), m_building.bottom_z - FT(10));
        const Point_3 p2 = Point_3(p.x(), p.y(), m_building.top_z + FT(10));
        line = Line_3(p1, p2);
        const FT z = internal::intersect_3(line, middle_plane);
        if (is_inside_building(z, p.z())) ++in;
        else ++out;
      }
    }

    bool is_inside_building(
      const FT curr_z, const FT real_z) const {

      return (
      curr_z > m_building.bottom_z - m_height_offset &&
      curr_z < real_z + m_height_offset );
    }
  };

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_INTERNAL_VISIBILITY_3_DEPR_H
