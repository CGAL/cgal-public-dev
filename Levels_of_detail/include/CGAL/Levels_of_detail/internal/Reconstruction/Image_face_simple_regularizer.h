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

#ifndef CGAL_LEVELS_OF_DETAIL_INTERNAL_IMAGE_FACE_SIMPLE_REGULARIZER_H
#define CGAL_LEVELS_OF_DETAIL_INTERNAL_IMAGE_FACE_SIMPLE_REGULARIZER_H

// STL includes.
#include <map>
#include <set>
#include <list>
#include <queue>
#include <vector>
#include <utility>
#include <memory>

// CGAL includes.
#include <CGAL/enum.h>
#include <CGAL/property_map.h>

// Internal includes.
#include <CGAL/Levels_of_detail/internal/utils.h>
#include <CGAL/Levels_of_detail/internal/struct.h>

// Other includes.
#include <CGAL/Levels_of_detail/internal/Reconstruction/Image.h>

// Testing.
#include "../../../../../test/Levels_of_detail/include/Saver.h"
#include "../../../../../test/Levels_of_detail/include/Utilities.h"

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

template<
typename GeomTraits,
typename Vertex,
typename Edge,
typename Halfedge,
typename Face,
typename Edge_type>
struct Image_face_simple_regularizer {

public:
  using Traits = GeomTraits;

  using FT = typename Traits::FT;
  using Point_2 = typename Traits::Point_2;
  using Point_3 = typename Traits::Point_3;
  using Segment_2 = typename Traits::Segment_2;
  using Segment_3 = typename Traits::Segment_3;
  using Line_2 = typename Traits::Line_2;
  using Line_3 = typename Traits::Line_3;
  using Vector_2 = typename Traits::Vector_2;
  using Plane_3 = typename Traits::Plane_3;
  using Triangle_2 = typename Traits::Triangle_2;
  using Triangle_3 = typename Traits::Triangle_3;
  using Intersect_2 = typename Traits::Intersect_2;
  using Intersect_3 = typename Traits::Intersect_3;

  using Image = internal::Image<Traits>;
  using Point_type = typename Image::Point_type;

  using FT_pair = std::pair<FT, FT>;
  using Indices = std::vector<std::size_t>;
  using Size_pair = std::pair<std::size_t, std::size_t>;

  using Saver = Saver<Traits>;
  using Inserter = Polygon_inserter<Traits>;
  using Indexer = internal::Indexer<Point_3>;
  using Color = CGAL::Color;

  Image_face_simple_regularizer(
    const std::vector<Segment_2>& boundary,
    std::vector<Vertex>& vertices,
    std::vector<Edge>& edges,
    std::vector<Halfedge>& halfedges,
    std::vector<Face>& faces,
    const Image& image,
    const std::map<std::size_t, Plane_3>& plane_map,
    const FT min_length_2,
    const FT angle_bound_2,
    const FT ordinate_bound_2) :
  m_boundary(boundary),
  m_vertices(vertices),
  m_edges(edges),
  m_halfedges(halfedges),
  m_faces(faces),
  m_image(image),
  m_plane_map(plane_map),
  m_min_length_2(min_length_2),
  m_angle_bound_2(angle_bound_2),
  m_ordinate_bound_2(ordinate_bound_2),
  m_pi(static_cast<FT>(CGAL_PI)),
  m_angle_threshold(FT(5)),
  m_bound_min(m_angle_bound_2),
  m_bound_max(FT(90) - m_bound_min)
  { }

  void set_face_edges(
    const std::vector<Edge>& segments) {

    clear();
    m_segments = segments;
  }

  void compute_multiple_directions(
    const Face& face) { }

  void regularize_face(
    Face& face) {

    CGAL_assertion(m_segments.size() != 0);
    /* if (m_angle_bound_2 == FT(0))
      return; */

    for (std::size_t i = 0; i < m_segments.size(); ++i) {
      const auto& edge = m_segments[i];
      const auto& from = m_vertices[edge.from_vertex];
      if (
        from.type == Point_type::CORNER ||
        from.type == Point_type::OUTER_CORNER ||
        from.type == Point_type::OUTER_BOUNDARY ) {
        simplify_adjacent_edges(i);
      }
    }
  }

  void snap(Face& face) {

    for (std::size_t i = 0; i < m_segments.size(); ++i) {
      const auto& edge = m_segments[i];
      const auto& from = m_vertices[edge.from_vertex];
      const auto& to   = m_vertices[edge.to_vertex];

      if (from.type == Point_type::OUTER_BOUNDARY) {
        snap_from(edge);
        continue;
      }

      if (to.type == Point_type::OUTER_BOUNDARY) {
        snap_to(edge);
        continue;
      }
    }
  }

  void clear() {
    m_segments.clear();
  }

private:
  const std::vector<Segment_2>& m_boundary;
  std::vector<Vertex>& m_vertices;
  std::vector<Edge>& m_edges;
  std::vector<Halfedge>& m_halfedges;
  std::vector<Face>& m_faces;
  const Image& m_image;
  const std::map<std::size_t, Plane_3>& m_plane_map;
  const FT m_min_length_2;
  const FT m_angle_bound_2;
  const FT m_ordinate_bound_2;
  const FT m_pi;
  const FT m_angle_threshold;
  const FT m_bound_min, m_bound_max;

  std::vector<Edge> m_segments;

  void snap_from(const Edge& edge) {

    const auto& segment = edge.segment;
    auto& from = m_vertices[edge.from_vertex];
    const auto& to = m_vertices[edge.to_vertex];
    const std::size_t bd_idx = from.bd_idx;
    const auto& bound = m_boundary[bd_idx];

    const FT angle   = angle_degree_2(segment, bound);
    const FT angle_2 = get_angle_2(angle);
    if (CGAL::abs(angle_2) >= m_bound_max) {
      const Line_2 line = Line_2(bound.source(), bound.target());
      auto proj = line.projection(to.point);

      if (CGAL::collinear_are_ordered_along_line(
        bound.source(), proj, bound.target())) {

        const Triangle_2 triangle = Triangle_2(from.point, to.point, proj);
        const FT area = CGAL::abs(triangle.area());
        if (area < FT(2))
          from.point = proj;
      }
    }

    const FT dist1 = internal::distance(bound.source(), from.point);
    const FT dist2 = internal::distance(from.point, bound.target());

    if (dist1 < dist2) {
      if (dist1 < m_ordinate_bound_2 * FT(2)) {
        from.point = bound.source();
      }
    } else {
      if (dist2 < m_ordinate_bound_2 * FT(2)) {
        from.point = bound.target();
      }
    }
  }

  void snap_to(const Edge& edge) {

    const auto& segment = edge.segment;
    const auto& from = m_vertices[edge.from_vertex];
    auto& to = m_vertices[edge.to_vertex];
    const std::size_t bd_idx = to.bd_idx;
    const auto& bound = m_boundary[bd_idx];

    const FT angle   = angle_degree_2(segment, bound);
    const FT angle_2 = get_angle_2(angle);
    if (CGAL::abs(angle_2) >= m_bound_max) {
      const Line_2 line = Line_2(bound.source(), bound.target());
      auto proj = line.projection(from.point);

      if (CGAL::collinear_are_ordered_along_line(
        bound.source(), proj, bound.target())) {

        const Triangle_2 triangle = Triangle_2(from.point, to.point, proj);
        const FT area = CGAL::abs(triangle.area());
        if (area < FT(2))
          to.point = proj;
      }
    }

    const FT dist1 = internal::distance(bound.source(), to.point);
    const FT dist2 = internal::distance(to.point, bound.target());

    if (dist1 < dist2) {
      if (dist1 < m_ordinate_bound_2 * FT(2)) {
        to.point = bound.source();
      }
    } else {
      if (dist2 < m_ordinate_bound_2 * FT(2)) {
        to.point = bound.target();
      }
    }
  }

  void simplify_adjacent_edges(
    const std::size_t start) {

    std::size_t curr = start;
    std::size_t idx  = m_segments[curr].from_vertex;
    Indices indices;
    apply_positive(start, curr, idx, indices);
    apply_negative(start, curr, idx, indices);
  }

  void apply_positive(
    const std::size_t start,
    std::size_t curr, std::size_t idx, Indices& indices) {

    Indices vs;
    indices.clear();
    const std::size_t n = m_segments.size();
    vs.push_back(m_segments[start].from_vertex);

    curr = (curr + 1) % n;
    idx = m_segments[curr].from_vertex;
    do {

      auto& vertex = m_vertices[idx];
      if (vertex.used) return;

      if (
        vertex.type == Point_type::CORNER ||
        vertex.type == Point_type::OUTER_CORNER ||
        vertex.type == Point_type::OUTER_BOUNDARY ) {

        vs.push_back(m_segments[curr].from_vertex);
        regularize_positive(start, indices, vs, curr);
        return;
      }

      indices.push_back(curr);
      vs.push_back(m_segments[curr].from_vertex);
      curr = (curr + 1) % n;
      idx = m_segments[curr].from_vertex;

    } while (curr != start);
  }

  void apply_negative(
    const std::size_t start,
    std::size_t curr, std::size_t idx, Indices& indices) {

    Indices vs;
    indices.clear();
    const std::size_t n = m_segments.size();
    vs.push_back(m_segments[start].from_vertex);

    curr = (curr + n - 1) % n;
    idx = m_segments[curr].from_vertex;
    do {

      auto& vertex = m_vertices[idx];
      if (vertex.used) return;

      if (
        vertex.type == Point_type::CORNER ||
        vertex.type == Point_type::OUTER_CORNER ||
        vertex.type == Point_type::OUTER_BOUNDARY ) {

        vs.push_back(m_segments[curr].from_vertex);
        regularize_negative(start, indices, vs, curr);
        return;
      }

      indices.push_back(curr);
      vs.push_back(m_segments[curr].from_vertex);
      curr = (curr + n - 1) % n;
      idx = m_segments[curr].from_vertex;

    } while (curr != start);
  }

  void regularize_positive(
    const std::size_t start,
    const Indices& indices,
    const Indices& vs,
    const std::size_t end) {

    if (indices.size() == 0)
      return;

    const FT eps = m_min_length_2 / FT(5);
    std::size_t curr = start;
    for (std::size_t i = 0; i < indices.size(); ++i) {

      const std::size_t next = indices[i];
      const auto& curr_edge = m_segments[curr];
      const auto& next_edge = m_segments[next];

      auto& vertex = m_vertices[next_edge.from_vertex];
      vertex.used = true;

      bool is_boundary = false;
      for (const std::size_t label : vertex.labels) {
        if (label == std::size_t(-1)) {
          is_boundary = true; break;
        }
      }
      if (is_boundary) continue;

      if (
        vertex.type == Point_type::LINEAR ||
        vertex.type == Point_type::FREE) {

        const auto p1 = get_prev_vertex(vs, vertex.index);
        const auto p2 = vertex.point;
        const auto p3 = get_next_vertex(vs, vertex.index);
        set_linear_point(p1, p2, p3, vertex);
      }
      curr = next;
    }
  }

  const Point_2 get_prev_vertex(
    const Indices& vs,
    const std::size_t query) {

    std::size_t idx = std::size_t(-1);
    for (std::size_t i = 0; i < vs.size(); ++i) {
      if (vs[i] == query) {
        idx = i; break;
      }
    }
    if (idx == std::size_t(-1)) {
      std::cout << "ERROR: error finding index prev!" << std::endl;
      return Point_2();
    }

    for (std::size_t i = idx - 1; i >= 0; i-=1) {
      if (!m_vertices[vs[i]].state)
        return m_vertices[vs[i]].point;
    }
    return m_vertices[vs[0]].point;
  }

  const Point_2 get_next_vertex(
    const Indices& vs,
    const std::size_t query) {

    std::size_t idx = std::size_t(-1);
    for (std::size_t i = 0; i < vs.size(); ++i) {
      if (vs[i] == query) {
        idx = i; break;
      }
    }
    if (idx == std::size_t(-1)) {
      std::cout << "ERROR: error finding index next!" << std::endl;
      return Point_2();
    }

    for (std::size_t i = idx + 1; i < vs.size(); i+=1) {
      if (!m_vertices[vs[i]].state)
        return m_vertices[vs[i]].point;
    }
    return m_vertices[vs[vs.size() - 1]].point;
  }

  void regularize_negative(
    const std::size_t start,
    const Indices& input,
    const Indices& vs,
    const std::size_t end) {

    if (input.size() == 0)
      return;

    auto indices = input;
    indices.push_back(end);

    const FT eps = m_min_length_2 / FT(5);
    for (std::size_t i = 0; i < indices.size() - 1; ++i) {
      const auto& curr_edge = m_segments[indices[i]];
      const auto& next_edge = m_segments[indices[i + 1]];

      auto& vertex = m_vertices[curr_edge.from_vertex];
      vertex.used = true;

      bool is_boundary = false;
      for (const std::size_t label : vertex.labels) {
        if (label == std::size_t(-1)) {
          is_boundary = true; break;
        }
      }
      if (is_boundary) continue;

      if (
        vertex.type == Point_type::LINEAR ||
        vertex.type == Point_type::FREE) {

        const auto p1 = get_prev_vertex(vs, vertex.index);
        const auto p2 = vertex.point;
        const auto p3 = get_next_vertex(vs, vertex.index);
        set_linear_point(p1, p2, p3, vertex);
      }
    }
  }

  void set_linear_point(
    const Point_2& p1, const Point_2& p2, const Point_2& p3,
    Vertex& vertex) {

    /*
    std::cout.precision(30);
    std::cout << p1 << std::endl;
    std::cout << p2 << std::endl;
    std::cout << p3 << std::endl; */

    const Triangle_2 triangle = Triangle_2(p1, p2, p3);
    const FT area = CGAL::abs(triangle.area());

    /* std::cout << area << std::endl << std::endl; */

    if (area < FT(2))
      vertex.state = true;
  }

  FT angle_degree_2(
    const Segment_2& longest, const Segment_2& segment) {

    const Vector_2 v1 =  segment.to_vector();
    const Vector_2 v2 = -longest.to_vector();

    const FT det = CGAL::determinant(v1, v2);
    const FT dot = CGAL::scalar_product(v1, v2);
    const FT angle_rad = static_cast<FT>(
      std::atan2(CGAL::to_double(det), CGAL::to_double(dot)));
    const FT angle_deg = angle_rad * FT(180) / m_pi;
    return angle_deg;
  }

  FT get_angle_2(const FT angle) {

    FT angle_2 = angle;
    if (angle_2 > FT(90)) angle_2 = FT(180) - angle_2;
    else if (angle_2 < -FT(90)) angle_2 = FT(180) + angle_2;
    return angle_2;
  }

  void save_face_contour(
    const std::vector<Edge>& edges,
    const std::size_t face_index) {

    std::vector<Segment_2> segments;
    segments.reserve(edges.size());
    for (const auto& edge : edges)
      segments.push_back(edge.segment);

    Saver saver;
    saver.save_polylines(
      segments, "/Users/monet/Documents/gf/lod/logs/buildings/tmp/contours/contour-" +
      std::to_string(face_index));
  }
};

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_INTERNAL_IMAGE_FACE_SIMPLE_REGULARIZER_H
