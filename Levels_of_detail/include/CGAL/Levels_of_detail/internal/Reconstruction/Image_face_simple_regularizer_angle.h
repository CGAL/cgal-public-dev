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

#ifndef CGAL_LEVELS_OF_DETAIL_INTERNAL_IMAGE_FACE_SIMPLE_REGULARIZER_ANGLE_H
#define CGAL_LEVELS_OF_DETAIL_INTERNAL_IMAGE_FACE_SIMPLE_REGULARIZER_ANGLE_H

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
struct Image_face_simple_regularizer_old {

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

  Image_face_simple_regularizer_angle(
    const std::vector<Segment_2>& boundary,
    std::vector<Vertex>& vertices,
    std::vector<Edge>& edges,
    std::vector<Halfedge>& halfedges,
    const Image& image,
    const std::map<std::size_t, Plane_3>& plane_map,
    const FT min_length_2,
    const FT angle_bound_2,
    const FT ordinate_bound_2) :
  m_boundary(boundary),
  m_vertices(vertices),
  m_edges(edges),
  m_halfedges(halfedges),
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

  void regularize_face(Face& face) {
    CGAL_assertion(m_segments.size() != 0);

    if (m_angle_bound_2 == FT(0))
      return;

    /* bool found = false; */
    for (std::size_t i = 0; i < m_segments.size(); ++i) {
      const auto& edge = m_segments[i];
      const auto& from = m_vertices[edge.from_vertex];
      if (
        from.type == Point_type::CORNER /* ||
        from.type == Point_type::OUTER_BOUNDARY ||
        from.type == Point_type::OUTER_CORNER */ ) {
        simplify_adjacent_edges(i);
      }

      /*
      if (
        from.type == Point_type::OUTER_BOUNDARY ||
        from.type == Point_type::OUTER_CORNER)
      found = true; */
    }

    /*
    if (!found) {
      face.skip = true;
      for (std::size_t i = 0; i < m_segments.size(); ++i) {
        const auto& edge = m_segments[i];
        auto& from = m_vertices[edge.from_vertex];
        from.state = true;
      }
    } */

    /* save_face_contour(m_segments, face.index); */
  }

  void clear() {
    m_segments.clear();
  }

private:
  const std::vector<Segment_2>& m_boundary;
  std::vector<Vertex>& m_vertices;
  std::vector<Edge>& m_edges;
  std::vector<Halfedge>& m_halfedges;
  const Image& m_image;
  const std::map<std::size_t, Plane_3>& m_plane_map;
  const FT m_min_length_2;
  const FT m_angle_bound_2;
  const FT m_ordinate_bound_2;
  const FT m_pi;
  const FT m_angle_threshold;
  const FT m_bound_min, m_bound_max;

  std::vector<Edge> m_segments;

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

    indices.clear();
    const std::size_t n = m_segments.size();

    curr = (curr + 1) % n;
    idx = m_segments[curr].from_vertex;
    do {

      auto& vertex = m_vertices[idx];
      if (vertex.used) return;

      if (
        vertex.type == Point_type::CORNER /* ||
        vertex.type == Point_type::OUTER_CORNER ||
        vertex.type == Point_type::OUTER_BOUNDARY */ ) {

        regularize_positive(start, indices);
        return;
      }

      indices.push_back(curr);
      curr = (curr + 1) % n;
      idx = m_segments[curr].from_vertex;

    } while (curr != start);
  }

  void apply_negative(
    const std::size_t start,
    std::size_t curr, std::size_t idx, Indices& indices) {

    indices.clear();
    const std::size_t n = m_segments.size();

    curr = (curr + n - 1) % n;
    idx = m_segments[curr].from_vertex;
    do {

      auto& vertex = m_vertices[idx];
      if (vertex.used) return;

      if (
        vertex.type == Point_type::CORNER /* ||
        vertex.type == Point_type::OUTER_CORNER ||
        vertex.type == Point_type::OUTER_BOUNDARY */ ) {

        const std::size_t end = (curr + n - 1) % n;
        regularize_negative(indices, end);
        return;
      }

      indices.push_back(curr);
      curr = (curr + n - 1) % n;
      idx = m_segments[curr].from_vertex;

    } while (curr != start);
  }

  void regularize_positive(
    const std::size_t start,
    const Indices& indices) {

    if (indices.size() == 0)
      return;

    const FT eps = m_min_length_2 / FT(5);
    std::size_t curr = start;
    for (const std::size_t next : indices) {

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

      const FT angle = angle_degree_2(curr_edge.segment, next_edge.segment);
      const FT angle_2 = get_angle_2(angle);
      if (CGAL::abs(angle_2) < m_bound_max) {
        vertex.state = true;
      } else {

        const FT dist1 = internal::distance(
          m_vertices[curr_edge.from_vertex].point,
          m_vertices[next_edge.from_vertex].point);
        const FT dist2 = internal::distance(
          m_vertices[next_edge.from_vertex].point,
          m_vertices[next_edge.to_vertex].point);
        if (dist1 < eps || dist2 < eps)
          vertex.state = true;
      }
      curr = next;
    }
  }

  void regularize_negative(
    const Indices& input,
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

      const FT angle = angle_degree_2(next_edge.segment, curr_edge.segment);
      const FT angle_2 = get_angle_2(angle);
      if (CGAL::abs(angle_2) < m_bound_max) {
        vertex.state = true;
      } else {

        const FT dist1 = internal::distance(
          m_vertices[next_edge.from_vertex].point,
          m_vertices[curr_edge.from_vertex].point);
        const FT dist2 = internal::distance(
          m_vertices[curr_edge.from_vertex].point,
          m_vertices[curr_edge.to_vertex].point);
        if (dist1 < eps || dist2 < eps)
          vertex.state = true;
      }
    }
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

#endif // CGAL_LEVELS_OF_DETAIL_INTERNAL_IMAGE_FACE_SIMPLE_REGULARIZER_ANGLE_H
