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

#ifndef CGAL_LEVELS_OF_DETAIL_INTERNAL_IMAGE_FACE_COMPLEX_REGULARIZER_H
#define CGAL_LEVELS_OF_DETAIL_INTERNAL_IMAGE_FACE_COMPLEX_REGULARIZER_H

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
struct Image_face_complex_regularizer {

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

  Image_face_complex_regularizer(
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
    const Face& face) {

    make_default_group(std::size_t(-1));
    get_multiple_directions();

    if (m_longest.size() != 0) {
      unify_along_contours();
      correct_directions();
    }
    if (m_longest.size() == 0)
      compute_longest_direction();

    /*
    std::cout << "Num face directions: " << m_longest.size() << std::endl;
    for (const std::size_t val : m_group)
      std::cout << val << " ";
    std::cout << std::endl; */
  }

  void regularize_face(const Face& face) {
    CGAL_assertion(m_segments.size() != 0);

    if (m_angle_bound_2 == FT(0))
      return;

    auto edges = m_segments;
    for (const auto& edge : edges) {
      auto& vertex = m_vertices[edge.from_vertex];
      if (!vertex.state) {
        vertex.used  = true;
        vertex.skip  = true;
      }
    }

    rotate_contour(edges);
    save_face_contour(edges, face.index);

    bool success = optimize_contour(edges);
    if (!success) {
      std::cout << "Error: face regularization failed, optimize!" << std::endl;
    } else {
      save_face_contour(edges, face.index);
    }

    success = connect_contour(edges);
    if (success) {
      save_face_contour(edges, face.index);

      for (auto& edge : edges) {
        auto& vertex = m_vertices[edge.from_vertex];
        if (!vertex.state) {
          vertex.skip  = false;
          vertex.point = edge.segment.source();
        }
      }
    } else {
      std::cout << "Error: face regularization failed, connect!" << std::endl;
    }

    for (const auto& edge : m_segments) {
      auto& vertex = m_vertices[edge.from_vertex];
      vertex.state = true;
    }
    m_segments = edges;
  }

  void clear() {
    m_segments.clear();
    m_bounds.clear();
    m_longest.clear();
    m_group.clear();
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
  std::vector<FT_pair> m_bounds;
  std::vector<Segment_2> m_longest;
  Indices m_group;

  void create_face_segments(
    const Face& face) {

    m_segments.clear();
    for (std::size_t i = 0; i < face.hedges.size(); ++i) {

      const std::size_t he_idx = face.hedges[i];
      const auto& he = m_halfedges[he_idx];
      const std::size_t from = he.from_vertex;
      const std::size_t to = he.to_vertex;

      const auto& s = m_vertices[from];
      const auto& t = m_vertices[to];

      if (!s.skip && !t.skip) {

        Edge edge;
        edge.from_vertex = from;
        edge.to_vertex = to;
        edge.segment = Segment_2(s.point, t.point);
        edge.type = m_edges[he.edg_idx].type;

        m_segments.push_back(edge);
        continue;
      }

      if (!s.skip && t.skip) {
        i = get_next(face, i);

        const std::size_t next_idx = face.hedges[i];
        const auto& next = m_halfedges[next_idx];
        const std::size_t end = next.to_vertex;
        const auto& other = m_vertices[end];

        Edge edge;
        edge.from_vertex = from;
        edge.to_vertex = end;
        edge.segment = Segment_2(s.point, other.point);
        edge.type = Edge_type::INTERNAL;

        m_segments.push_back(edge);
        continue;
      }
    }
  }

  std::size_t get_next(
    const Face& face,
    const std::size_t start) {

    for (std::size_t i = start; i < face.hedges.size(); ++i) {
      const std::size_t he_idx = face.hedges[i];
      const auto& he = m_halfedges[he_idx];
      const std::size_t to = he.to_vertex;
      if (m_vertices[to].skip) continue;
      return i;
    }

    const std::size_t i = face.hedges.size() - 1;
    const std::size_t he_idx = face.hedges[i];
    const auto& he = m_halfedges[he_idx];
    const std::size_t to = he.to_vertex;
    return i;
  }

  void make_default_group(
    const std::size_t value) {

    CGAL_assertion(m_segments.size() != 0);

    m_group.clear();
    m_group.resize(m_segments.size(), value);
  }

  void get_multiple_directions() {

    Segment_2 linear_segment;
    const bool line_found = find_linear_segment(linear_segment);

    if (line_found) {
      m_longest.push_back(linear_segment);
      m_bounds.push_back(std::make_pair(FT(45), FT(45)));
    }

    std::vector<Segment_2> boundaries;
    const bool boundary_found = find_boundary_segments(boundaries);

    if (boundary_found) {
      for (const auto& segment : boundaries) {
        m_longest.push_back(segment);
        m_bounds.push_back(std::make_pair(FT(45), FT(45)));
      }
    }

    if (m_longest.size() == 0) return;
    create_group();
  }

  bool find_linear_segment(
    Segment_2& linear_segment) {

    for (const auto& edge : m_segments) {
      if (edge.type == Edge_type::BOUNDARY)
        continue;

      const auto& from = m_vertices[edge.from_vertex];
      const auto& to = m_vertices[edge.to_vertex];

      if (
        from.type == Point_type::LINEAR &&
        to.type == Point_type::LINEAR) {
        linear_segment = edge.segment;
        return true;
      }
    }
    return false;
  }

  bool find_boundary_segments(
    std::vector<Segment_2>& boundaries) {

    boundaries.clear();
    for (const auto& edge : m_segments) {
      if (edge.type == Edge_type::BOUNDARY) {
        boundaries.push_back(edge.segment);
        continue;
      }
    }

    std::sort(boundaries.begin(), boundaries.end(),
    [](const Segment_2& a, const Segment_2& b) {
      return a.squared_length() > b.squared_length();
    });

    return ( boundaries.size() != 0 );
  }

  void create_group() {

    for (std::size_t i = 0; i < m_segments.size(); ++i) {
      auto& edge = m_segments[i];
      edge.skip = false;

      const auto& from = m_vertices[edge.from_vertex];
      const auto& to = m_vertices[edge.to_vertex];

      /*
      if (
        from.skip || to.skip) {

        std::cout << "should not be here" << std::endl;
        edge.skip = true;
        m_group[i] = std::size_t(-1); continue;
      } */

      if (from.used && to.used)
        edge.skip = true;

      if (from.type == Point_type::LINEAR && to.type == Point_type::LINEAR)
        edge.skip = true;

      if (edge.type == Edge_type::BOUNDARY)
        edge.skip = true;

      for (std::size_t j = 0; j < m_longest.size(); ++j) {
        const FT angle = angle_degree_2(m_longest[j], edge.segment);
        const FT angle_2 = get_angle_2(angle);

        if (
          (CGAL::abs(angle_2) <= m_bound_min) ||
          (CGAL::abs(angle_2) >= m_bound_max) )  {

          m_group[i] = j; break;
        }
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

  void unify_along_contours() {

    for (std::size_t i = 0; i < m_segments.size(); ++i) {
      const auto& edge = m_segments[i];
      if (edge.skip) continue;

      if (m_group[i] == std::size_t(-1)) {

        const std::size_t m = m_segments.size();
        std::size_t im = (i + m - 1) % m;
        std::size_t ip = (i + 1) % m;

        bool stop = false;
        std::size_t max_count = 0;
        do {

          if (m_group[im] != std::size_t(-1)) {
            m_group[i] = m_group[im]; break;
          }

          if (m_group[ip] != std::size_t(-1)) {
            m_group[i] = m_group[ip]; break;
          }

          im = (im + m - 1) % m;
          ip = (ip + 1) % m;

          if (im == i || ip == i) stop = true;
          ++max_count;

        } while (!stop && max_count < m * 2);
        if (stop || max_count > m * 2)
          m_group[i] = 0;
      }
    }
  }

  void correct_directions() {

    const std::size_t n = m_segments.size();
    Indices group; group.reserve(n);
    for (std::size_t i = 0; i < n; ++i) {

      const std::size_t im = (i + n - 1) % n;
      const std::size_t ip = (i + 1) % n;

      const std::size_t gm = m_group[im];
      const std::size_t gi = m_group[i];
      const std::size_t gp = m_group[ip];

      if (gm == gp && gi != gm)
        group.push_back(gm);
      else
        group.push_back(gi);
    }
    m_group = group;
  }

  void compute_longest_direction() {

    m_bounds.clear();
    m_longest.clear();
    m_group.clear();
    make_default_group(std::size_t(-1));

    std::vector<Segment_2> segments;
    segments.reserve(m_segments.size());

    for (auto& edge : m_segments) {
      edge.skip = false;
      segments.push_back(edge.segment);
    }

    std::sort(segments.begin(), segments.end(),
    [](const Segment_2& a, const Segment_2& b) {
      return a.squared_length() > b.squared_length();
    });

    CGAL_assertion(segments.size() != 0);
    m_longest.push_back(segments[0]);
    m_bounds.push_back(std::make_pair(FT(45), FT(45)));

    for (std::size_t i = 0; i < m_segments.size(); ++i) {
      auto& edge = m_segments[i];

      const auto& from = m_vertices[edge.from_vertex];
      const auto& to = m_vertices[edge.to_vertex];

      if (from.used && to.used)
        edge.skip = true;

      /*
      if (
        from.skip || to.skip) {
        edge.skip = true;
        m_group[i] = std::size_t(-1);
      } else {
        edge.skip = false;
        m_group[i] = 0;
      } */

      m_group[i] = 0;
    }
  }

  void rotate_contour(std::vector<Edge>& edges) {

    for (std::size_t i = 0; i < edges.size(); ++i) {
      auto& edge = edges[i];
      if (edge.skip) continue;

      const std::size_t gr_idx = m_group[i];
      if (gr_idx == std::size_t(-1))
        continue;

      auto& segment = edge.segment;
      const auto& longest_segment = m_longest[gr_idx];
      const auto& bounds = m_bounds[gr_idx];

      const bool success = rotate_segment(longest_segment, bounds, segment);
      if (!success)
        m_group[i] = std::size_t(-1);
    }
  }

  bool rotate_segment(
    const Segment_2& longest_segment,
    const FT_pair& bounds,
    Segment_2& segment) {

    const FT angle = angle_degree_2(longest_segment, segment);
    const FT angle_2 = get_angle_2(angle);

    if (CGAL::abs(angle_2) <= bounds.first)
      rotate(angle, FT(180), longest_segment, segment); // parallel case
    if (CGAL::abs(angle_2) >= bounds.second)
      rotate(angle, FT(90), longest_segment, segment); // orthogonal case
    return true;
  }

  void rotate(
    const FT angle_2,
    const FT ref_angle_2,
    const Segment_2& longest_segment,
    Segment_2& segment) {

    FT angle = angle_2;
    if (angle < FT(0)) angle = angle + ref_angle_2;
    else if (angle > FT(0)) angle = angle - ref_angle_2;

    Point_2 source = segment.source();
    Point_2 target = segment.target();

    const Point_2 b = internal::middle_point_2(source, target);
    const FT angle_rad = angle * m_pi / FT(180);

    internal::rotate_point_2(angle_rad, b, source);
    internal::rotate_point_2(angle_rad, b, target);

    segment = Segment_2(source, target);
  }

  bool optimize_contour(std::vector<Edge>& edges) {

    std::vector<Edge> clean;
    remove_zero_length_segments(edges, clean);
    if (clean.size() < 4)
      return false;

    std::vector< std::vector<Edge> > groups;
    create_consecutive_groups(clean, groups);

    /*
    std::size_t num_groups = 0;
    for (const auto& group : groups)
      if (group.size() > 1) ++num_groups;
    std::cout << "Num consecutive groups: " << num_groups << std::endl; */

    for (auto& group : groups)
      if (group.size() > 1)
        optimize_group(group);

    edges.clear();
    for (const auto& group : groups) {
      for (const auto& edge : group)
        edges.push_back(edge);
    }
    return true;
  }

  void remove_zero_length_segments(
    const std::vector<Edge>& edges,
    std::vector<Edge>& clean) {

    clean.clear();
    for (const auto& edge : edges)
      if (edge.segment.squared_length() > internal::tolerance<FT>())
        clean.push_back(edge);
  }

  void create_consecutive_groups(
    const std::vector<Edge>& edges,
    std::vector< std::vector<Edge> >& groups) {

    groups.clear();
    const std::size_t num_edges = edges.size();
    std::vector<bool> states(num_edges, false);

    std::vector<Edge> group;
    std::size_t gr_idx = 0;

    for (std::size_t i = 0; i < num_edges; ++i) {
      const auto& edge_i = edges[i];
      if (states[i]) continue;

      group.clear();
      group.push_back(edge_i);

      if (edge_i.skip) {
        groups.push_back(group);
        ++gr_idx; continue;
      }
      states[i] = true;

      const std::size_t ip = (i + 1) % num_edges;
      if (ip != 0) {
        for (std::size_t j = ip; j < num_edges; ++j) {
          const auto& edge_j = edges[j];
          if (edge_j.skip) break;

          const FT angle = angle_degree_2(
            edge_i.segment, edge_j.segment);
          const FT angle_2 = get_angle_2(angle);

          if (CGAL::abs(angle_2) <= m_angle_threshold) {
            group.push_back(edge_j); states[j] = true;
          } else break;
        }
      }
      groups.push_back(group);
      ++gr_idx;
    }
  }

  void optimize_group(
    std::vector<Edge>& edges) {

    Indices indices;
    indices.reserve(edges.size());
    for (std::size_t i = 0; i < edges.size(); ++i)
      indices.push_back(i);

    // Sort.
    std::sort(indices.begin(), indices.end(),
    [&edges](const std::size_t idx1, const std::size_t idx2) -> bool {
      const FT length_1 = edges[idx1].segment.squared_length();
      const FT length_2 = edges[idx2].segment.squared_length();
      return length_1 > length_2;
    });

    std::vector< std::vector<Edge> > groups;
    std::map<std::size_t, std::size_t> seg_map;
    std::vector<bool> states(edges.size(), false);

    std::vector<Edge> group;
    std::size_t gr_idx = 0;

    // Create groups.
    for (std::size_t i = 0; i < indices.size(); ++i) {
      const int idxi = static_cast<int>(indices[i]);
      const auto& edge_i = edges[idxi];
      if (states[idxi]) continue;

      group.clear();
      group.push_back(edge_i);
      seg_map[idxi] = gr_idx;
      states[idxi] = true;

      const auto p =
      internal::middle_point_2(
        edge_i.segment.source(), edge_i.segment.target());

      // Go right.
      const int idxip = idxi + 1;
      if (idxi < edges.size() - 1 && !states[idxip]) {
        for (int j = idxip; j < edges.size(); ++j) {
          const auto& edge_j = edges[j];
          if (states[j]) break;

          const Line_2 line = Line_2(
            edge_j.segment.source(), edge_j.segment.target());
          const auto q = line.projection(p);
          const FT distance = internal::distance(p, q);

          if (distance <= m_ordinate_bound_2) {
            group.push_back(edge_j); states[j] = true;
            seg_map[j] = gr_idx;
          } else break;
        }
      }

      // Go left.
      const int idxim = idxi - 1;
      if (idxi > 0 && !states[idxim]) {
        auto j = idxim;
        while (j >= 0) {
          const auto& edge_j = edges[j];
          if (states[j]) break;

          const Line_2 line = Line_2(
            edge_j.segment.source(), edge_j.segment.target());
          const auto q = line.projection(p);
          const FT distance = internal::distance(p, q);

          if (distance <= m_ordinate_bound_2) {
            group.push_back(edge_j); states[j] = true;
            seg_map[j] = gr_idx;
          } else break;
          --j;
        }
      }

      groups.push_back(group);
      ++gr_idx;
    }

    // Create lines.
    std::vector<Line_2> lines;
    lines.reserve(groups.size());
    for (const auto& group : groups) {

      const Segment_2 segment = find_weighted_segment(group);
      const Line_2 line = Line_2(segment.source(), segment.target());
      lines.push_back(line);
    }

    // Project.
    Point_2 p, q;
    for (std::size_t i = 0; i < edges.size(); ++i) {
      const std::size_t gr_idx = seg_map.at(i);
      const Line_2& line = lines[gr_idx];

      auto& edge = edges[i];

      const auto& s = edge.segment.source();
      const auto& t = edge.segment.target();

      p = line.projection(s);
      q = line.projection(t);
      edge.segment = Segment_2(p, q);
    }

    // Final result.
    std::vector<Edge> result;
    for (std::size_t i = 0; i < edges.size(); ++i) {
      const std::size_t ip = (i + 1) % edges.size();
      if (ip == 0) {
        result.push_back(edges[i]);
        break;
      }

      const auto gri = seg_map.at(i);
      const auto grj = seg_map.at(ip);

      result.push_back(edges[i]);
      if (gri != grj) {

        const auto& edgei = edges[i];
        const auto& edgej = edges[ip];

        Line_2 line = Line_2(
          edgej.segment.source(), edgej.segment.target());
        auto source = internal::middle_point_2(
          edgei.segment.source(), edgei.segment.target());
        auto target = line.projection(source);
        Segment_2 orth = Segment_2(source, target);

        Edge extra;
        extra.from_vertex = edgei.to_vertex;
        extra.to_vertex = edgej.from_vertex;
        extra.segment = orth;
        extra.type = Edge_type::INTERNAL;
        extra.skip = false;

        result.push_back(extra);
      }
    }
    edges = result;
  }

  bool connect_contour(std::vector<Edge>& edges) {

    bool success = false;

    success = clean_segments(edges);
    if (!success) return false;

    success = make_segments_collinear(edges);
    if (!success) return false;

    intersect_segments(edges);

    /*
    success = clean_and_intersect_segments(edges);
    if (success)
      success = clean_segments(edges); */
    return success;
  }

  bool clean_segments(
    std::vector<Edge>& edges) {

    std::vector<Edge> clean;
    remove_zero_length_segments(edges, clean);
    if (clean.size() < 4)
      return false;

    std::vector<Edge> filtered;
    const bool success = filter_out_wrong_segments(clean, filtered);
    if (filtered.size() < 4 || !success)
      return false;

    edges = filtered;
    return true;
  }

  bool filter_out_wrong_segments(
    const std::vector<Edge>& edges,
    std::vector<Edge>& filtered) {

    filtered.clear();
    const std::size_t n = edges.size();
    const std::size_t start = find_initial_index(edges);

    std::size_t i = start;
    std::vector<Edge> parallel_edges;
    std::size_t max_count = 0;
    do {

      const std::size_t prev = i;
      const bool success = get_parallel_segments(
        edges, parallel_edges, i);
      CGAL_assertion(parallel_edges.size() != 0);
      if (!success) return false;
      const std::size_t curr = i;

      Segment_2 segment;
      const FT sum_length =
      create_segment_from_parallel_segments(parallel_edges, segment);

      const auto& prev_edge = edges[prev];
      const auto& curr_edge = edges[curr];

      if (parallel_edges.size() > 1) {

        Edge extra;
        extra.from_vertex = prev_edge.to_vertex;
        extra.to_vertex = curr_edge.from_vertex;
        extra.segment = segment;
        extra.type = Edge_type::INTERNAL;
        extra.skip = false;
        filtered.push_back(extra);

      } else if (parallel_edges.size() == 1) {
        filtered.push_back(parallel_edges[0]);
      }

      ++max_count;
    } while (i != start && max_count < n * 2);
    if (max_count > n * 2) return false;
    return true;
  }

  std::size_t find_initial_index(
    const std::vector<Edge>& edges) {

    const std::size_t n = edges.size();
    for (std::size_t i = 0; i < n; ++i) {
      if (edges[i].skip) continue;

      const std::size_t im = (i + n - 1) % n;
      const std::size_t ip = (i + 1) % n;

      const auto& si = edges[i].segment;
      const auto& sm = edges[im].segment;
      const auto& sp = edges[ip].segment;

      const auto pair = is_parallel_segment(sm, si, sp);
      const bool previous_is_orthogonal = !(pair.first);
      if (previous_is_orthogonal) return i;
    }
    return 0;
  }

  std::pair<bool, bool> is_parallel_segment(
    const Segment_2& sm, const Segment_2& si, const Segment_2& sp) {

    const FT angle_mi = angle_degree_2(sm, si);
    const FT angle_mi_2 = get_angle_2(angle_mi);

    const FT angle_pi = angle_degree_2(si, sp);
    const FT angle_pi_2 = get_angle_2(angle_pi);

    const bool source_cond = ( CGAL::abs(angle_mi_2) <= m_angle_threshold );
    const bool target_cond = ( CGAL::abs(angle_pi_2) <= m_angle_threshold );

    return std::make_pair(source_cond, target_cond);
  }

  bool get_parallel_segments(
    const std::vector<Edge>& edges,
    std::vector<Edge>& parallel_edges,
    std::size_t& seed) {

    parallel_edges.clear();
    const std::size_t n = edges.size();

    std::size_t i = seed;
    bool next_is_parallel = false;
    std::size_t max_count = 0;
    do {

      const std::size_t im = (i + n - 1) % n;
      const std::size_t ip = (i + 1) % n;

      const auto& si = edges[i].segment;
      const auto& sm = edges[im].segment;
      const auto& sp = edges[ip].segment;

      parallel_edges.push_back(edges[i]);
      if (edges[i].skip || edges[ip].skip) {
        seed = ip; return true;
      }

      const auto pair = is_parallel_segment(sm, si, sp);
      next_is_parallel = pair.second;
      i = ip;

      ++max_count;
    } while (next_is_parallel && max_count < n * 2);
    if (max_count > n * 2) return false;
    seed = i;
    return true;
  }

  FT create_segment_from_parallel_segments(
    const std::vector<Edge>& parallel_edges,
    Segment_2& result) {

    Segment_2 ref_segment = find_weighted_segment(parallel_edges);
    const Line_2 line =
    Line_2(ref_segment.source(), ref_segment.target());

    FT sum_length = FT(0);
    std::vector<Point_2> points;
    for (const auto& edge : parallel_edges) {

      const Point_2 p = line.projection(edge.segment.source());
      const Point_2 q = line.projection(edge.segment.target());

      points.push_back(p);
      points.push_back(q);

      sum_length += internal::distance(p, q);
    }
    update_segment(points, ref_segment);
    result = ref_segment;
    return sum_length;
  }

  Segment_2 find_weighted_segment(
    const std::vector<Edge>& edges) {

    std::vector<FT> weights;
    compute_distance_weights(edges, weights);
    const Segment_2 ref_segment = find_central_segment(edges);
    const Segment_2 result =
      compute_weighted_segment(edges, weights, ref_segment);

    if (result.source() == result.target())
      return ref_segment;
    return result;
  }

  void compute_distance_weights(
    const std::vector<Edge>& edges,
    std::vector<FT>& weights) {

    weights.clear();
    weights.reserve(edges.size());

    FT sum_distance = FT(0);
    for (const auto& edge : edges) {
      const FT distance =
        internal::distance(
          edge.segment.source(), edge.segment.target());
      sum_distance += distance;
      weights.push_back(distance);
    }

    for (auto& weight : weights)
      weight /= sum_distance;
  }

  Segment_2 find_central_segment(
    const std::vector<Edge>& edges) {

    std::vector<Edge> tmp;
    for (const auto& edge : edges) {
      if (edge.type == Edge_type::BOUNDARY)
        tmp.push_back(edge);
    }
    if (tmp.size() != 0)
      return find_longest_segment(tmp);

    tmp.clear();
    for (const auto& edge : edges) {
      if (edge.skip)
        tmp.push_back(edge);
    }
    if (tmp.size() != 0)
      return find_longest_segment(tmp);

    Point_2 source, target;
    FT x1 = FT(0), y1 = FT(0);
    FT x2 = FT(0), y2 = FT(0);
    for (const auto& edge : edges) {
      x1 += edge.segment.source().x();
      x2 += edge.segment.target().x();

      y1 += edge.segment.source().y();
      y2 += edge.segment.target().y();
    }

    const FT size = static_cast<FT>(edges.size());
    x1 /= size; y1 /= size;
    x2 /= size; y2 /= size;

    source = Point_2(x1, y1);
    target = Point_2(x2, y2);

    if (source == target)
      return find_longest_segment(edges);
    return Segment_2(source, target);
  }

  Segment_2 find_longest_segment(
    const std::vector<Edge>& edges) {

    FT max_length = -FT(1);
    std::size_t edg_idx = std::size_t(-1);

    for (std::size_t i = 0; i < edges.size(); ++i) {
      const auto& edge = edges[i];
      const FT length = edge.segment.squared_length();
      if (length > max_length) {

        max_length = length;
        edg_idx = i;
      }
    }
    return edges[edg_idx].segment;
  }

  Segment_2 compute_weighted_segment(
    const std::vector<Edge>& edges,
    const std::vector<FT>& weights,
    const Segment_2& ref_segment) {

    const Point_2& s = ref_segment.source();
    const Point_2& t = ref_segment.target();

    const Point_2 b = internal::middle_point_2(s, t);

    Vector_2 dir = Vector_2(FT(0), FT(0));
    for (std::size_t i = 0; i < weights.size(); ++i) {
      const FT weight = weights[i];

      const Segment_2& segment = edges[i].segment;
      const Line_2 line = Line_2(segment.source(), segment.target());
      const Point_2 p = line.projection(b);

      const Vector_2 v = Vector_2(b, p);
      dir += v * weight;
    }

    const Point_2 news = s + dir;
    const Point_2 newt = t + dir;

    return Segment_2(news, newt);
  }

  void update_segment(
    const std::vector<Point_2>& points,
    Segment_2& segment) {

    FT min_proj_value =  internal::max_value<FT>();
    FT max_proj_value = -internal::max_value<FT>();

    const Vector_2 ref_vector = segment.to_vector();
    Point_2 ref_point;
    internal::compute_barycenter_2(points, ref_point);

    Point_2 p, q;
    for (const auto& point : points) {
      const Vector_2 curr_vector(ref_point, point);
      const FT value = CGAL::scalar_product(curr_vector, ref_vector);

      if (value < min_proj_value) {
        min_proj_value = value;
        p = point; }
      if (value > max_proj_value) {
        max_proj_value = value;
        q = point; }
    }
    segment = Segment_2(p, q);
  }

  bool make_segments_collinear(
    std::vector<Edge>& edges) {

    std::map<std::size_t, std::size_t> seg_map;
    std::vector< std::vector<Edge> > groups;
    create_collinear_groups(edges, groups, seg_map);

    /*
    std::size_t num_groups = 0;
    for (const auto& group : groups)
      if (group.size() > 1) ++num_groups;
    std::cout << "Num collinear groups: " << num_groups << std::endl; */

    std::vector<Line_2> lines;
    lines.reserve(groups.size());
    for (const auto& group : groups) {

      const Segment_2 segment = find_central_segment(group);
      const Line_2 line = Line_2(segment.source(), segment.target());
      lines.push_back(line);
    }

    for (std::size_t i = 0; i < edges.size(); ++i) {
      auto& edge = edges[i];
      if (edge.skip) continue; // comment it out if you want to project all edges

      const std::size_t gr_idx = seg_map.at(i);
      const Line_2& line = lines[gr_idx];

      const auto& s = edge.segment.source();
      const auto& t = edge.segment.target();

      // Do not use return here. It will not work! Tested!
      if (line.a() == FT(0) && line.b() == FT(0) && line.c() == FT(0))
        continue;

      const Point_2 p = line.projection(s);
      const Point_2 q = line.projection(t);

      edge.segment = Segment_2(p, q);
    }
    return true;
  }

  void create_collinear_groups(
    const std::vector<Edge>& edges,
    std::vector< std::vector<Edge> >& groups,
    std::map<std::size_t, std::size_t>& seg_map) {

    groups.clear(); seg_map.clear();
    std::vector<bool> states(edges.size(), false);

    std::vector<Edge> group;
    std::size_t gr_idx = 0;

    for (std::size_t i = 0; i < edges.size(); ++i) {
      const auto& edge_i = edges[i];
      if (states[i]) continue;

      group.clear();
      group.push_back(edge_i);
      seg_map[i] = gr_idx;
      states[i] = true;

      const auto p = internal::middle_point_2(
        edge_i.segment.source(), edge_i.segment.target());
      for (std::size_t j = 0; j < edges.size(); ++j) {
        const auto& edge_j = edges[j];
        if (states[j]) continue;

        const FT angle   = angle_degree_2(edge_i.segment, edge_j.segment);
        const FT angle_2 = get_angle_2(angle);

        if (CGAL::abs(angle_2) <= m_angle_threshold) {
          const Line_2 line = Line_2(
            edge_j.segment.source(), edge_j.segment.target());

          const auto q = line.projection(p);
          const FT distance = internal::distance(p, q);

          if (distance <= m_ordinate_bound_2) {
            group.push_back(edge_j);
            seg_map[j] = gr_idx;
            states[j] = true;
          }
        }
      }
      groups.push_back(group);
      ++gr_idx;
    }
  }

  void intersect_segments(
    std::vector<Edge>& edges) {

    std::vector<Edge> intersected;
    const std::size_t n = edges.size();
    for (std::size_t i = 0; i < n; ++i) {
      const std::size_t ip = (i + 1) % n;

      auto& edi = edges[i];
      auto& edp = edges[ip];

      auto& si = edi.segment;
      auto& sp = edp.segment;

      if (!edi.skip) {
        if (!edp.skip) {
          intersect_segment(si, sp);
          intersected.push_back(edi);
        } else {
          intersected.push_back(edi);
          if (!internal::are_equal_points_2(si.target(), sp.source())) {
            Edge extra;
            extra.from_vertex = edi.to_vertex;
            extra.to_vertex = edp.from_vertex;
            extra.segment = Segment_2(si.target(), sp.source());
            extra.type = Edge_type::INTERNAL;
            extra.skip = false;
            intersected.push_back(extra);
          }
        }
      } else {
        intersected.push_back(edi);
        if (!internal::are_equal_points_2(si.target(), sp.source())) {
          Edge extra;
          extra.from_vertex = edi.to_vertex;
          extra.to_vertex = edp.from_vertex;
          extra.segment = Segment_2(si.target(), sp.source());
          extra.type = Edge_type::INTERNAL;
          extra.skip = false;
          intersected.push_back(extra);
        }
      }
    }
    edges = intersected;
  }

  void intersect_segment(
    Segment_2& si, Segment_2& sp) {

    const Line_2 line_1 = Line_2(si.source(), si.target());
    const Line_2 line_2 = Line_2(sp.source(), sp.target());

    Point_2 source, target, point;
    const bool success = intersect_2(line_1, line_2, point);
    if (success) {

      source = si.source();
      target = si.target();
      target = point;
      si = Segment_2(source, target);

      source = sp.source();
      target = sp.target();
      source = point;
      sp = Segment_2(source, target);
    }
  }

  bool intersect_2(
    const Line_2& line_1, const Line_2& line_2,
    Point_2& in_point) {

    typename std::result_of<Intersect_2(Line_2, Line_2)>::type result
    = CGAL::intersection(line_1, line_2);
    if (result) {
      if (const Line_2* line = boost::get<Line_2>(&*result))
        return false;
      else {
        const Point_2* point = boost::get<Point_2>(&*result);
        in_point = *point; return true;
      }
    }
    return false;
  }

  bool clean_and_intersect_segments(
    std::vector<Edge>& edges) {

    const bool success = clean_segments(edges);
    if (!success) return false;

    intersect_segments(edges);

    for (const auto& edge : edges) {
      const auto& s = edge.segment.source();
      const auto& t = edge.segment.target();

      if (
        std::isnan(CGAL::to_double(s.x())) || std::isnan(CGAL::to_double(s.y())) ||
        std::isnan(CGAL::to_double(t.x())) || std::isnan(CGAL::to_double(t.y())) )
      return false;
    }
    return true;
  }

  void save_face_contour(
    const std::vector<Edge>& edges,
    const std::size_t face_index) {

    std::vector<Segment_2> segments;
    segments.reserve(edges.size());

    for (const auto& edge : edges) {
      segments.push_back(edge.segment);

      /* segments.push_back(Segment_2(
        m_vertices[edge.from_vertex].point, m_vertices[edge.to_vertex].point)); */
    }

    Saver saver;
    saver.save_polylines(
      segments, "/Users/monet/Documents/gf/lod/logs/buildings/tmp/contours/contour-" +
      std::to_string(face_index));
  }
};

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_INTERNAL_IMAGE_FACE_COMPLEX_REGULARIZER_H
