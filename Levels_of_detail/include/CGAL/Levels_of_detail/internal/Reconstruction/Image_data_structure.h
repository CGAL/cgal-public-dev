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

#ifndef CGAL_LEVELS_OF_DETAIL_INTERNAL_IMAGE_DATA_STRUCTURE_H
#define CGAL_LEVELS_OF_DETAIL_INTERNAL_IMAGE_DATA_STRUCTURE_H

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
#include <CGAL/point_generators_3.h>
#include <CGAL/Random.h>

// Internal includes.
#include <CGAL/Levels_of_detail/internal/utils.h>
#include <CGAL/Levels_of_detail/internal/struct.h>

// Other includes.
#include <CGAL/Levels_of_detail/internal/Reconstruction/Image.h>
#include <CGAL/Levels_of_detail/internal/Shape_detection/Region_growing.h>
#include <CGAL/Levels_of_detail/internal/Spatial_search/K_neighbor_query.h>
#include <CGAL/Levels_of_detail/internal/Reconstruction/Image_face_simplifier.h>
#include <CGAL/Levels_of_detail/internal/Reconstruction/Image_face_simple_regularizer.h>
#include <CGAL/Levels_of_detail/internal/Reconstruction/Image_face_complex_regularizer.h>

// Testing.
#include "../../../../../test/Levels_of_detail/include/Saver.h"
#include "../../../../../test/Levels_of_detail/include/Utilities.h"

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

template<typename GeomTraits>
struct Image_data_structure {

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
  using Vector_3 = typename Traits::Vector_3;
  using Plane_3 = typename Traits::Plane_3;
  using Triangle_2 = typename Traits::Triangle_2;
  using Triangle_3 = typename Traits::Triangle_3;
  using Intersect_2 = typename Traits::Intersect_2;
  using Intersect_3 = typename Traits::Intersect_3;

  using Image = internal::Image<Traits>;

  using My_point = typename Image::My_point;
  using Contour = typename Image::Contour;
  using Point_type = typename Image::Point_type;
  using Pixels = typename Image::Pixels;
  using Pixel_point_map = typename Image::Pixel_point_map;

  using Indices = std::vector<std::size_t>;
  using Size_pair = std::pair<std::size_t, std::size_t>;
  using Triangle_set_3 = std::vector<Triangle_3>;

  using Saver = Saver<Traits>;
  using Inserter = Polygon_inserter<Traits>;
  using Indexer = internal::Indexer<Point_3>;
  using Color = CGAL::Color;

  using K_neighbor_query =
  internal::K_neighbor_query<Traits, Pixels, Pixel_point_map>;
  using Triangulation =
  internal::Triangulation<Traits>;
  using LF_circulator =
  typename Triangulation::Delaunay::Line_face_circulator;

  using Random = CGAL::Random;

  struct Vertex {
    Point_2 point;
    std::size_t index = std::size_t(-1);
    std::size_t bd_idx = std::size_t(-1);
    Point_type type = Point_type::DEFAULT;
    std::set<std::size_t> labels;
    Indices hedges;
    Indices neighbors;
    Indices faces;

    bool used  = false;
    bool skip  = false;
    bool state = false;

    std::set<std::size_t> tmp;

    void clear() {
      labels.clear();
      point = Point_2();
      index = std::size_t(-1);
      bd_idx = std::size_t(-1);
      type = Point_type::DEFAULT;
      hedges.clear();
      neighbors.clear();
    }
  };

  enum class Edge_type {
    DEFAULT  = 0,
    BOUNDARY = 1,
    INTERNAL = 2
  };

  struct Edge {
    std::size_t index = std::size_t(-1);
    std::size_t from_vertex = std::size_t(-1);
    std::size_t to_vertex = std::size_t(-1);
    Size_pair labels = std::make_pair(std::size_t(-1), std::size_t(-1));
    Segment_2 segment;
    Edge_type type = Edge_type::DEFAULT;
    Size_pair faces = std::make_pair(std::size_t(-1), std::size_t(-1));
    FT length = FT(0);
    double weight = 1.0;

    bool used = false;
    bool skip = false;

    double get_length() const {
      return static_cast<double>(length);
    }

    FT compute_length() const {
      return internal::distance(segment.source(), segment.target());
    }
  };

  struct Halfedge {
    std::size_t index = std::size_t(-1);
    std::size_t edg_idx = std::size_t(-1);
    std::size_t next = std::size_t(-1);
    std::size_t opposite = std::size_t(-1);
    std::size_t from_vertex = std::size_t(-1);
    std::size_t to_vertex = std::size_t(-1);
  };

  enum class Face_type {
    DEFAULT = 0,
    CLOSED = 1,
    BOUNDARY = 2,
    INTERNAL = 3
  };

  struct Face {
    std::size_t index = std::size_t(-1);
    Face_type type = Face_type::DEFAULT;
    std::size_t label = std::size_t(-1);
    std::set<std::size_t> probs; // label probabilities
    Indices hedges;
    Indices neighbors;
    Triangulation tri;
    FT area = FT(0);
    double weight = 1.0;
    std::size_t original = std::size_t(-1);

    bool used = false;
    bool skip = false;

    std::size_t level = std::size_t(-1);
    std::set<std::size_t> tmp;

    double get_area() const {
      return static_cast<double>(area);
    }
  };

  using Edges = std::vector<Edge>;

  class DS_neighbor_query {

  public:
    DS_neighbor_query(
      const std::vector<Vertex>& vertices,
      const std::vector<Edge>& edges,
      const std::vector<Halfedge>& halfedges) :
    m_vertices(vertices),
    m_edges(edges),
    m_halfedges(halfedges)
    { }

    void operator()(
      const std::size_t query_index,
      Indices& neighbors) const {

      neighbors.clear();
      const auto& he = m_halfedges[query_index];
      if (he.next != std::size_t(-1))
        neighbors.push_back(he.next);
    }

  private:
    const std::vector<Vertex>& m_vertices;
    const std::vector<Edge>& m_edges;
    const std::vector<Halfedge>& m_halfedges;
  };

  class DS_region_type {

  public:
    DS_region_type(
      const std::vector<Vertex>& vertices,
      const std::vector<Edge>& edges,
      const std::vector<Halfedge>& halfedges) :
    m_vertices(vertices),
    m_edges(edges),
    m_halfedges(halfedges)
    { }

    bool is_already_visited(
      const std::size_t, const std::size_t, const bool) const {
      return false;
    }

    bool is_part_of_region(
      const std::size_t, const std::size_t,
      const Indices&) const {
      return true;
    }

    bool is_valid_region(const Indices& region) const {
      return region.size() >= 3;
    }

    void update(const Indices&) {
      // skipped!
    }

  private:
    const std::vector<Vertex>& m_vertices;
    const std::vector<Edge>& m_edges;
    const std::vector<Halfedge>& m_halfedges;
  };

  using Region_growing = internal::Region_growing<
    std::vector<Halfedge>, DS_neighbor_query, DS_region_type>;

  Image_data_structure(
    std::vector<Segment_2>& boundary,
    const std::vector<Image>& ridges,
    const Image& image,
    const std::map<std::size_t, Plane_3>& plane_map,
    const FT noise_level_2,
    const FT min_length_2,
    const FT angle_bound_2,
    const FT ordinate_bound_2,
    const FT max_height_difference,
    const FT top_z) :
  m_boundary(boundary),
  m_ridges(ridges),
  m_image(image),
  m_plane_map(plane_map),
  m_noise_level_2(noise_level_2),
  m_min_length_2(min_length_2),
  m_angle_bound_2(angle_bound_2),
  m_ordinate_bound_2(ordinate_bound_2),
  m_max_height_difference(max_height_difference),
  m_top_z(top_z),
  m_pi(static_cast<FT>(CGAL_PI)),
  m_random(0),
  m_knq(
    m_image.pixels,
    FT(24),
    m_image.pixel_point_map)
  { }

  std::vector<Vertex>& vertices() {
    return m_vertices;
  }

  std::vector<Edge>& edges() {
    return m_edges;
  }

  std::vector<Halfedge>& halfedges() {
    return m_halfedges;
  }

  std::vector<Face>& faces() {
    return m_faces;
  }

  void build() {

    initialize_vertices();
    initialize_edges();
    initialize_faces();

    default_vertex_states();
    snap_to_boundary();
    for (auto& face : m_faces) {
      const bool success = update_face(face);
      if (!success) face.skip = true;
    }
    default_vertex_states();

    CGAL_assertion(
      m_halfedges.size() == m_edges.size() * 2);

    /*
    std::cout.precision(30);
    std::cout << "num vertices: " << m_vertices.size() << std::endl;
    for (const auto& vertex : m_vertices) {
      if (
        vertex.type != Point_type::FREE &&
        vertex.type != Point_type::LINEAR) {

        if (vertex.type != Point_type::CORNER)
          continue;

        std::cout <<
        int(vertex.type) << " : " <<
        vertex.bd_idx << " , " <<
        vertex.labels.size() << " , " <<
        vertex.neighbors.size() << " , " <<
        vertex.hedges.size() << std::endl;
      }
    } */

    /*
    std::cout << "num edges: " << m_edges.size() << std::endl;
    std::cout << "num halfedges: " << m_halfedges.size() << std::endl;
    for (const auto& edge : m_edges) {
      std::cout <<
        edge.labels.first << " " <<
        edge.labels.second << std::endl;
    } */

    /*
    std::cout << "num edges: " << m_edges.size() << std::endl;
    std::cout << "num halfedges: " << m_halfedges.size() << std::endl;
    for (const auto& he : m_halfedges) {
      std::cout <<
      he.from_vertex << " " << he.to_vertex << " : " <<
      int(m_vertices[he.from_vertex].type) << " " <<
      int(m_vertices[he.to_vertex].type) << " : " <<
      he.edg_idx << " " <<
      he.index << " " << he.opposite << std::endl;
    } */

    /*
    std::cout << "num faces: " << m_faces.size() << std::endl;
    for (const auto& face : m_faces) {
      std::cout <<
      int(face.type) << " : " <<
      face.label << " , " <<
      face.probs.size() << " , " <<
      face.neighbors.size() << " , " <<
      face.hedges.size() << std::endl;
    } */

    /*
    save_faces_polylines("initial");
    save_faces_ply("initial"); */

    std::cout << "data structure built" << std::endl;
  }

  void simplify() {

    default_vertex_states();
    using Image_face_simplifier = internal::Image_face_simplifier<
      Traits, Vertex, Edge, Halfedge, Face>;
    Image_face_simplifier image_face_simplifier(
      m_boundary,
      m_vertices, m_edges, m_halfedges,
      m_image, m_plane_map,
      m_noise_level_2);

    for (const auto& face : m_faces) {
      if (face.skip) continue;
      image_face_simplifier.simplify_face(face);
    }

    make_skip();
    default_vertex_states();
    for (auto& face : m_faces) {
      const bool success = update_face(face);
      if (!success) face.skip = true;
    }

    /*
    save_faces_polylines("simplified");
    save_faces_ply("simplified"); */

    save_all_faces_ply(1000, "regularized");
  }

  void regularize(const std::size_t level) {
    regularize_simple(level);
  }

  void regularize_simple(const std::size_t level) {

    default_vertex_states();
    using Image_face_regularizer = internal::Image_face_simple_regularizer<
      Traits, Vertex, Edge, Halfedge, Face, Edge_type>;
    Image_face_regularizer image_face_regularizer(
      m_boundary,
      m_vertices, m_edges, m_halfedges, m_faces,
      m_image, m_plane_map,
      m_min_length_2, m_angle_bound_2, m_ordinate_bound_2);

    std::vector<Edges> face_edges(m_faces.size());
    for (std::size_t i = 0; i < m_faces.size(); ++i) {
      const auto& face = m_faces[i];
      auto& edges = face_edges[i];
      create_face_edges(face, edges);
    }
    update_edge_neighbors(face_edges);

    for (std::size_t i = 0; i < m_faces.size(); ++i) {
      auto& face = m_faces[i];
      if (face.skip) continue;
      image_face_regularizer.set_face_edges(face_edges[i]);
      image_face_regularizer.compute_multiple_directions(face);
      image_face_regularizer.regularize_face(face);
    }

    make_skip();
    default_vertex_states();
    for (auto& face : m_faces) {
      const bool success = update_face(face);
      if (!success) face.skip = true;
    }
    mark_bad_faces();

    /*
    save_faces_polylines("faces");
    save_faces_ply("faces"); */

    save_all_faces_ply(level, "regularized");
  }

  void project_linear(const std::size_t level) {

    std::vector<Edges> face_edges(m_faces.size());
    for (std::size_t i = 0; i < m_faces.size(); ++i) {
      const auto& face = m_faces[i];
      auto& edges = face_edges[i];
      create_face_edges(face, edges);
    }
    update_edge_neighbors(face_edges);

    for (auto& vertex : m_vertices) {
      if (vertex.skip) continue;
      project_onto_line(face_edges, vertex);
    }

    for (auto& face : m_faces) {
      const bool success = update_face(face);
      if (!success) face.skip = true;
    }
    save_all_faces_ply(level, "regularized");
  }

  void project_onto_line(
    const std::vector<Edges>& face_edges,
    Vertex& vertex) {

    Indices labels;
    find_vertex_labels(face_edges, vertex, labels);
    if (labels.size() < 2) return;

    const std::size_t l0 = labels[0];
    const std::size_t l1 = labels[1];

    const auto& plane0 = m_plane_map.at(l0);
    const auto& plane1 = m_plane_map.at(l1);

    Line_2 line;
    bool success = intersect_planes(plane0, plane1, line);
    if (!success) return;

    if (vertex.type == Point_type::LINEAR) {
      const auto proj = line.projection(vertex.point);
      if (is_close_by(vertex, proj))
        vertex.point = proj;
      return;
    }

    if (vertex.type == Point_type::OUTER_BOUNDARY) {
      if (vertex.bd_idx == std::size_t(-1)) return;

      const auto& segment = m_boundary[vertex.bd_idx];
      const auto other = Line_2(segment.source(), segment.target());

      Point_2 res;
      success = intersect_2(line, other, res);
      if (!success) return;

      if (is_close_by(vertex, res)) {
        if (CGAL::collinear_are_ordered_along_line(
            segment.source(), res, segment.target())) {
          vertex.point = res;
        } else {

          const FT dist1 = internal::distance(res, segment.source());
          const FT dist2 = internal::distance(res, segment.target());
          if (dist1 <= dist2) {
            vertex.point = segment.source();
          } else {
            vertex.point = segment.target();
          }
        }
      }
      return;
    }

    if (vertex.type == Point_type::CORNER) {
      if (labels.size() == 3) {
        const auto& plane1 = m_plane_map.at(labels[0]);
        const auto& plane2 = m_plane_map.at(labels[1]);
        const auto& plane3 = m_plane_map.at(labels[2]);

        Point_2 res;
        const bool success = intersect_planes(plane1, plane2, plane3, res);
        if (!success) return;

        if (is_close_by(vertex, res)) {
          if (is_inside_building(labels, res))
            vertex.point = res;
        }
        return;
      }

      std::vector<Point_2> points;
      for (std::size_t i = 0; i < labels.size(); ++i) {
        for (std::size_t j = 0; j < labels.size(); ++j) {
          if (j <= i) continue;

          const std::size_t li = labels[i];
          const std::size_t lj = labels[j];

          const auto& planei = m_plane_map.at(li);
          const auto& planej = m_plane_map.at(lj);

          Line_2 other;
          bool success = intersect_planes(planei, planej, other);
          if (!success) continue;

          const auto proj = other.projection(vertex.point);
          if (is_close_by(vertex, proj)) {
            if (is_inside_building(labels, proj))
              points.push_back(proj);
          }
        }
      }

      FT min_dist = internal::max_value<FT>();
      std::size_t idx = std::size_t(-1);
      for (std::size_t i = 0; i < points.size(); ++i) {
        const FT dist = internal::distance(points[i], vertex.point);
        if (dist < min_dist) {
          min_dist = dist; idx = i;
        }
      }
      if (idx != std::size_t(-1))
        vertex.point = points[idx];
      return;
    }
  }

  bool intersect_planes(
    const Plane_3& plane1,
    const Plane_3& plane2,
    const Plane_3& plane3,
    Point_2& point_2) {

    auto result = CGAL::intersection(plane1, plane2, plane3);
    Point_3 point_3; bool found = false;
    if (result) {
      if (const Point_3* p = boost::get<Point_3>(&*result)) {
        found = true; point_3 = *p;
      }
    }
    if (!found) return false;

    point_2 = Point_2(point_3.x(), point_3.y());
    return true;
  }

  bool is_close_by(
    const Vertex& vertex,
    const Point_2& query) {

    const FT distance = internal::distance(vertex.point, query);
    return distance < m_noise_level_2 * FT(3);
  }

  bool is_inside_building(
    const Indices& labels, const Point_2& query) {

    for (const auto& face : m_faces) {
      if (face.skip) continue;
      bool found = false;
      for (const std::size_t label : labels) {
        if (face.label == label) {
          found = true; break;
        }
      }
      if (!found) continue;

      const auto& delaunay = face.tri.delaunay;
      const auto fh = delaunay.locate(query);
      if (fh->info().tagged)
        return true;
    }
    return false;
  }

  void find_vertex_labels(
    const std::vector<Edges>& face_edges,
    const Vertex& query,
    Indices& labels) {

    std::set<std::size_t> labs;
    for (std::size_t i = 0; i < face_edges.size(); ++i) {
      for (std::size_t j = 0; j < face_edges[i].size(); ++j) {
        const auto& edge = face_edges[i][j];
        if (
          edge.from_vertex == query.index ||
          edge.to_vertex == query.index) {

          const std::size_t f1 = edge.faces.first;
          const std::size_t f2 = edge.faces.second;

          if (f1 != std::size_t(-1)) {
            const std::size_t l1 = m_faces[f1].label;
            if (l1 != std::size_t(-1))
              labs.insert(l1);
          }

          if (f2 != std::size_t(-1)) {
            const std::size_t l2 = m_faces[f2].label;
            if (l2 != std::size_t(-1))
              labs.insert(l2);
          }
        }
      }
    }

    labels.clear();
    for (const std::size_t lab : labs)
      labels.push_back(lab);
  }

  void snap(const std::size_t level) {

    using Image_face_regularizer = internal::Image_face_simple_regularizer<
      Traits, Vertex, Edge, Halfedge, Face, Edge_type>;
    Image_face_regularizer image_face_regularizer(
      m_boundary,
      m_vertices, m_edges, m_halfedges, m_faces,
      m_image, m_plane_map,
      m_min_length_2, m_angle_bound_2, m_ordinate_bound_2);

    std::vector<Edges> face_edges(m_faces.size());
    for (std::size_t i = 0; i < m_faces.size(); ++i) {
      const auto& face = m_faces[i];
      auto& edges = face_edges[i];
      create_face_edges(face, edges);
    }
    update_edge_neighbors(face_edges);

    for (std::size_t i = 0; i < m_faces.size(); ++i) {
      auto& face = m_faces[i];
      if (face.skip) continue;
      image_face_regularizer.set_face_edges(face_edges[i]);
      image_face_regularizer.snap(face);
    }

    make_skip();

    default_vertex_states();
    for (auto& face : m_faces) {
      const bool success = update_face(face);
      if (!success) face.skip = true;
    }

    mark_bad_faces();
    save_all_faces_ply(level, "regularized");
  }

  void merge_corners(const std::size_t level) {

    std::vector<Edges> face_edges(m_faces.size());
    for (std::size_t i = 0; i < m_faces.size(); ++i) {
      const auto& face = m_faces[i];
      auto& edges = face_edges[i];
      create_face_edges(face, edges);
    }
    update_edge_neighbors(face_edges);

    std::set<std::size_t> faces;
    std::set<std::size_t> vs;
    Indices vec;

    for (auto& vertex : m_vertices) {
      if (vertex.type == Point_type::CORNER) {

        faces.clear();
        for (const std::size_t he_idx : vertex.hedges) {
          const auto& edge = m_edges[m_halfedges[he_idx].edg_idx];

          if (edge.faces.first != std::size_t(-1))
            faces.insert(edge.faces.first);
          if (edge.faces.second != std::size_t(-1))
            faces.insert(edge.faces.second);
        }

        vs.clear();
        for (const std::size_t fi : faces)
          add_vertices(vertex, face_edges[fi], vs);

        vec.clear();
        for (const std::size_t vidx : vs)
          vec.push_back(vidx);

        std::sort(vec.begin(), vec.end(),
        [&](const std::size_t i1, const std::size_t i2) -> bool {
          const FT dist1 = internal::distance(vertex.point, m_vertices[i1].point);
          const FT dist2 = internal::distance(vertex.point, m_vertices[i2].point);

          return dist1 < dist2;
        });

        if (vec.size() == 3) {

          auto& p = vertex.point;
          const auto& q1 = m_vertices[vec[0]].point;
          const auto& q2 = m_vertices[vec[1]].point;

          const Triangle_2 triangle = Triangle_2(p, q1, q2);
          const FT area = CGAL::abs(triangle.area());
          if (area < FT(3)) {
            const Line_2 line1 = Line_2(p, m_vertices[vec[2]].point);
            const Line_2 line2 = Line_2(q1, q2);
            Point_2 res;
            intersect_2(line1, line2, res);

            Indices labels;
            find_vertex_labels(face_edges, vertex, labels);
            if (is_inside_building(labels, res))
              p = res;
          }
        }
      }
    }

    for (auto& face : m_faces) {
      const bool success = update_face(face);
      if (!success) face.skip = true;
    }
    save_all_faces_ply(level, "regularized");
  }

  void merge_free_parts(const std::size_t level) {

    std::vector<Edges> face_edges(m_faces.size());
    for (std::size_t i = 0; i < m_faces.size(); ++i) {
      const auto& face = m_faces[i];
      auto& edges = face_edges[i];
      create_face_edges(face, edges);
    }
    update_edge_neighbors(face_edges);

    std::set<std::size_t> faces;
    std::set<std::size_t> vs;
    Indices vec;

    for (auto& vertex : m_vertices) {
      if (
        vertex.type == Point_type::FREE ||
        vertex.type == Point_type::LINEAR) {

        faces.clear();
        for (const std::size_t he_idx : vertex.hedges) {
          const auto& edge = m_edges[m_halfedges[he_idx].edg_idx];

          if (edge.faces.first != std::size_t(-1))
            faces.insert(edge.faces.first);
          if (edge.faces.second != std::size_t(-1))
            faces.insert(edge.faces.second);
        }

        vs.clear();
        for (const std::size_t fi : faces)
          add_vertices(vertex, face_edges[fi], vs);

        vec.clear();
        for (const std::size_t vidx : vs)
          vec.push_back(vidx);

        if (vec.size() == 2) {

          std::size_t type1 = std::size_t(-1);
          std::size_t bd_idx1 = std::size_t(-1);
          const FT angle1 = get_smallest_angle(
            vertex, m_vertices[vec[0]], type1, bd_idx1);

          std::size_t type2 = std::size_t(-1);
          std::size_t bd_idx2 = std::size_t(-1);
          const FT angle2 = get_smallest_angle(
            vertex, m_vertices[vec[1]], type2, bd_idx2);

          if (bd_idx1 == std::size_t(-1) && bd_idx2 == std::size_t(-1))
            continue;

          if (bd_idx1 != std::size_t(-1) && bd_idx2 == std::size_t(-1)) {
            if (type1 == 0) { // parallel
              const Line_2 line = Line_2(
                m_boundary[bd_idx1].source(), m_boundary[bd_idx1].target());
              vertex.point = line.projection(vertex.point);
            } else { // orthogonal
              const Line_2 line = Line_2(
                m_boundary[bd_idx1].source(), m_boundary[bd_idx1].target());
              const auto orth = line.perpendicular(m_vertices[vec[0]].point);
              vertex.point = orth.projection(vertex.point);
            }
            continue;
          }

          if (bd_idx1 == std::size_t(-1) && bd_idx2 != std::size_t(-1)) {
            if (type2 == 0) { // parallel
              const Line_2 line = Line_2(
                m_boundary[bd_idx2].source(), m_boundary[bd_idx2].target());
              vertex.point = line.projection(vertex.point);
            } else { // orthogonal
              const Line_2 line = Line_2(
                m_boundary[bd_idx2].source(), m_boundary[bd_idx2].target());
              const auto orth = line.perpendicular(m_vertices[vec[1]].point);
              vertex.point = orth.projection(vertex.point);
            }
            continue;
          }

          FT diff1, diff2;
          if (type1 == 0) diff1 = FT(45) - CGAL::abs(get_angle_2(angle1));
          else diff1 = FT(90) - CGAL::abs(get_angle_2(angle1));
          if (type2 == 0) diff2 = FT(45) - CGAL::abs(get_angle_2(angle2));
          else diff2 = FT(90) - CGAL::abs(get_angle_2(angle2));

          if (diff1 < diff2) {
            if (type1 == 0) { // parallel
              const Line_2 line = Line_2(
                m_boundary[bd_idx1].source(), m_boundary[bd_idx1].target());
              vertex.point = line.projection(vertex.point);
            } else { // orthogonal
              const Line_2 line = Line_2(
                m_boundary[bd_idx1].source(), m_boundary[bd_idx1].target());
              const auto orth = line.perpendicular(m_vertices[vec[0]].point);
              vertex.point = orth.projection(vertex.point);
            }
          } else {
            if (type2 == 0) { // parallel
              const Line_2 line = Line_2(
                m_boundary[bd_idx2].source(), m_boundary[bd_idx2].target());
              vertex.point = line.projection(vertex.point);
            } else { // orthogonal
              const Line_2 line = Line_2(
                m_boundary[bd_idx2].source(), m_boundary[bd_idx2].target());
              const auto orth = line.perpendicular(m_vertices[vec[1]].point);
              vertex.point = orth.projection(vertex.point);
            }
          }
        }
      }
    }

    for (auto& face : m_faces) {
      const bool success = update_face(face);
      if (!success) face.skip = true;
    }
    save_all_faces_ply(level, "regularized");
  }

  FT get_smallest_angle(
    const Vertex& query, const Vertex& other,
    std::size_t& type, std::size_t& bd_idx) {

    if (other.type == Point_type::OUTER_BOUNDARY) {
      bd_idx = other.bd_idx;

      const auto segment = Segment_2(query.point, other.point);
      const FT angle = angle_degree_2(
        m_boundary[bd_idx], segment);
      const FT angle_2 = CGAL::abs(get_angle_2(angle));
      if (angle_2 <= FT(45)) type = 0;
      else type = 1;
      return angle;
    }

    if (other.type == Point_type::OUTER_CORNER) {

      Indices bdis;
      const std::size_t n = m_boundary.size();
      bdis.push_back((other.bd_idx + n - 1) % n);
      bdis.push_back(other.bd_idx);

      const auto segment = Segment_2(query.point, other.point);
      FT min_diff = FT(1000000);

      FT a = 0;
      for (std::size_t i = 0; i < bdis.size(); ++i) {
        const std::size_t bdi = bdis[i];

        const FT angle = angle_degree_2(m_boundary[bdi], segment);
        const FT angle_2 = CGAL::abs(get_angle_2(angle));

        FT diff = 0;
        std::size_t t = std::size_t(-1);

        if (angle_2 <= FT(45)) {
          diff = FT(45) - angle_2;
          t = 0;
        } else {
          diff = FT(90) - angle_2;
          t = 1;
        }

        if (diff <= min_diff) {
          min_diff = diff;
          type = t;
          bd_idx = bdi;
          a = angle;
        }
      }
      return a;
    }
    return FT(0);
  }

  void create_all_face_edges() {

    m_face_edges.clear();
    m_face_edges.resize(m_faces.size());

    for (std::size_t i = 0; i < m_faces.size(); ++i) {
      const auto& face = m_faces[i];
      auto& edges = m_face_edges[i];
      create_face_edges(face, edges);
    }
    update_edge_neighbors(m_face_edges);
  }

  void get_wall_outer_segments(
    std::vector<Segment_3>& segments_3) {

    segments_3.clear();
    for (std::size_t i = 0; i < m_faces.size(); ++i) {
      const auto& face = m_faces[i];
      if (face.skip) continue;
      const auto& edges = m_face_edges[i];

      Plane_3 plane;
      if (face.label != std::size_t(-1))
        plane = m_plane_map.at(face.label);
      else
        plane = Plane_3(
          Point_3(FT(0), FT(0), m_top_z),
          Vector_3(FT(0), FT(0), FT(1)));

      for (const auto& edge : edges) {
        if (edge.type == Edge_type::BOUNDARY) {
          const auto& s = edge.segment.source();
          const auto& t = edge.segment.target();

          /*
          const Point_3 p = internal::position_on_plane_3(s, plane);
          const Point_3 q = internal::position_on_plane_3(t, plane); */

          const Point_3 p = get_position_on_plane_3(edge.from_vertex, s, plane);
          const Point_3 q = get_position_on_plane_3(edge.to_vertex, t, plane);

          segments_3.push_back(Segment_3(p, q));
        }
      }
    }
  }

  void get_wall_inner_segments(
    std::vector<Segment_3>& segments_3) {

    segments_3.clear();
    for (std::size_t i = 0; i < m_faces.size(); ++i) {
      const auto& face = m_faces[i];
      if (face.skip) continue;
      const auto& edges = m_face_edges[i];

      Plane_3 plane;
      if (face.label != std::size_t(-1))
        plane = m_plane_map.at(face.label);
      else
        plane = Plane_3(
          Point_3(FT(0), FT(0), m_top_z),
          Vector_3(FT(0), FT(0), FT(1)));

      for (const auto& edge : edges) {
        if (edge.type == Edge_type::INTERNAL) {

          const std::size_t f1 = edge.faces.first;
          const std::size_t f2 = edge.faces.second;

          if (f1 != std::size_t(-1) && f2 != std::size_t(-1)) {
            const auto& face1 = m_faces[f1];
            const auto& face2 = m_faces[f2];

            if (face1.label == face2.label)
              continue;
          }

          const auto& s = edge.segment.source();
          const auto& t = edge.segment.target();

          /*
          const Point_3 p = internal::position_on_plane_3(s, plane);
          const Point_3 q = internal::position_on_plane_3(t, plane); */

          const Point_3 p = get_position_on_plane_3(edge.from_vertex, s, plane);
          const Point_3 q = get_position_on_plane_3(edge.to_vertex, t, plane);

          segments_3.push_back(Segment_3(p, q));
          add_points(edge, f1, f2);
        }
      }
    }
  }

  void add_points(
    const Edge& edge,
    const std::size_t f1, const std::size_t f2) {

    return;
    const auto& face1 = m_faces[f1];
    const auto& face2 = m_faces[f2];

    const auto& segment = edge.segment;
    const auto& s = segment.source();
    const auto& t = segment.target();

    const auto& plane1 = m_plane_map.at(face1.label);
    const auto& plane2 = m_plane_map.at(face2.label);

    const auto p1 = internal::position_on_plane_3(s, plane1);
    const auto p2 = internal::position_on_plane_3(t, plane1);

    const auto q1 = internal::position_on_plane_3(s, plane2);
    const auto q2 = internal::position_on_plane_3(t, plane2);

    const Triangle_3 triangle1 = Triangle_3(p1, p2, q2);
    const Triangle_3 triangle2 = Triangle_3(q2, q1, p1);

    if (CGAL::sqrt(triangle1.squared_area()) > FT(1) / FT(10))
      add_points(triangle1);
    if (CGAL::sqrt(triangle2.squared_area()) > FT(1) / FT(10))
      add_points(triangle2);
  }

  void add_points(
    const Triangle_3& triangle) {

    std::ofstream outfile;
    outfile.precision(20);
    outfile.open(
      "/Users/monet/Documents/lod/logs/buildings/02_building_interior_points.ply",
      std::ios_base::app);

    std::vector<Point_3> samples;
    using Point_generator = CGAL::Random_points_in_triangle_3<Point_3>;

    Point_generator generator(triangle, m_random);
    std::copy_n(
      generator, 200, std::back_inserter(samples));
    for (const auto& sample : samples)
      outfile << sample << " 255 102 51 " << std::endl;
    outfile.close();
  }

  void get_roof_triangles(
    std::vector<Triangle_set_3>& triangle_sets_3) {

    triangle_sets_3.clear();
    std::set<std::size_t> labels;
    for (const auto& face : m_faces) {
      if (face.skip) continue;
      labels.insert(face.label);
    }

    Triangulation tri;
    Triangle_set_3 triangle_set_3;
    for (const std::size_t label : labels) {

      Plane_3 plane;
      if (label != std::size_t(-1))
        plane = m_plane_map.at(label);
      else
        plane = Plane_3(
          Point_3(FT(0), FT(0), m_top_z),
          Vector_3(FT(0), FT(0), FT(1)));

      tri.delaunay.clear();
      for (std::size_t i = 0; i < m_faces.size(); ++i) {
        const auto& face = m_faces[i];
        if (face.skip) continue;
        if (face.label != label) continue;

        const auto& edges = m_face_edges[i];
        for (const auto& edge : edges) {

          const std::size_t f1 = edge.faces.first;
          const std::size_t f2 = edge.faces.second;

          if (f1 != std::size_t(-1) && f2 != std::size_t(-1)) {
            const auto& face1 = m_faces[f1];
            const auto& face2 = m_faces[f2];

            if (face1.label == face2.label)
              continue;
          }

          const auto& s = edge.segment.source();
          const auto& t = edge.segment.target();

          const auto vh1 = tri.delaunay.insert(s);
          const auto vh2 = tri.delaunay.insert(t);

          vh1->info().object_index = edge.from_vertex;
          vh2->info().object_index = edge.to_vertex;

          if (vh1 != vh2)
            tri.delaunay.insert_constraint(vh1, vh2);
        }
      }

      Face tmp_face;
      tmp_face.tri = tri;
      create_face_visibility(tmp_face);

      triangle_set_3.clear();
      for (auto fh = tmp_face.tri.delaunay.finite_faces_begin();
      fh != tmp_face.tri.delaunay.finite_faces_end(); ++fh) {
        if (!fh->info().interior) continue;

        const auto& p0 = fh->vertex(0)->point();
        const auto& p1 = fh->vertex(1)->point();
        const auto& p2 = fh->vertex(2)->point();

        const Point_3 q0 = get_position_on_plane_3(
          fh->vertex(0)->info().object_index, p0, plane);
        const Point_3 q1 = get_position_on_plane_3(
          fh->vertex(1)->info().object_index, p1, plane);
        const Point_3 q2 = get_position_on_plane_3(
          fh->vertex(2)->info().object_index, p2, plane);

        triangle_set_3.push_back(Triangle_3(q0, q1, q2));
      }
      triangle_sets_3.push_back(triangle_set_3);
    }
  }

  void save_faces_polylines(const std::string folder) {

    std::vector<Edge> edges;
    std::vector<Segment_2> segments;

    for (const auto& face : m_faces) {
      if (face.skip) continue;

      create_face_edges(face, edges);

      segments.clear();
      for (const auto& edge : edges) {
        segments.push_back(edge.segment);
      }

      Saver saver;
      saver.save_polylines(
        segments, "/Users/monet/Documents/gf/lod/logs/buildings/tmp/" +
        folder + "/image-faces-" +
        std::to_string(face.index));
    }
  }

  void save_faces_ply(const std::string folder) {

    for (const auto& face : m_faces) {
      if (face.skip) continue;

      const FT z = FT(0);
      std::size_t num_vertices = 0;
      Indexer indexer;

      std::vector<Point_3> vertices;
      std::vector<Indices> faces;
      std::vector<Color> fcolors;

      Inserter inserter(faces, fcolors);
      auto output_vertices = std::back_inserter(vertices);
      auto output_faces = boost::make_function_output_iterator(inserter);
      face.tri.output_with_label_color(
        indexer, num_vertices, output_vertices, output_faces, z);

      Saver saver;
      saver.export_polygon_soup(
        vertices, faces, fcolors,
        "/Users/monet/Documents/gf/lod/logs/buildings/tmp/" +
        folder + "/faces-" +
        std::to_string(face.index));
    }
  }

  void save_all_faces_ply(
    const std::size_t level,
    const std::string folder) {

    const FT z = FT(0);
    std::size_t num_vertices = 0;
    Indexer indexer;

    std::vector<Point_3> vertices;
    std::vector<Indices> faces;
    std::vector<Color> fcolors;

    Inserter inserter(faces, fcolors);
    auto output_vertices = std::back_inserter(vertices);
    auto output_faces = boost::make_function_output_iterator(inserter);

    for (const auto& face : m_faces) {
      if (face.skip) continue;
      face.tri.output_with_label_color(
        indexer, num_vertices, output_vertices, output_faces, z);
    }

    Saver saver;
    saver.export_polygon_soup(
      vertices, faces, fcolors,
      "/Users/monet/Documents/gf/lod/logs/buildings/tmp/" +
      folder + "/tree-level-" +
      std::to_string(level));
  }

  void clear() {
    m_vertices.clear();
    m_edges.clear();
    m_faces.clear();
    m_halfedges.clear();
  }

private:
  std::vector<Segment_2>& m_boundary;
  const std::vector<Image>& m_ridges;
  const Image& m_image;
  const std::map<std::size_t, Plane_3>& m_plane_map;

  const FT m_noise_level_2;
  const FT m_min_length_2;
  const FT m_angle_bound_2;
  const FT m_ordinate_bound_2;
  const FT m_max_height_difference;
  const FT m_top_z;
  const FT m_pi;

  Random m_random;
  K_neighbor_query m_knq;

  std::vector<Vertex>   m_vertices;
  std::vector<Edge>     m_edges;
  std::vector<Halfedge> m_halfedges;
  std::vector<Face>     m_faces;

  std::map<Point_2, std::size_t> m_vertex_map;
  std::map<Size_pair, Halfedge> m_halfedge_map;
  std::map<std::size_t, std::size_t> m_face_map;
  std::vector<Edges> m_face_edges;
  std::map<size_t, Size_pair> m_ends;

  void add_vertices(
    const Vertex& query,
    const Edges& edges,
    std::set<std::size_t>& vs) {

    std::size_t idx = std::size_t(-1);
    const std::size_t n = edges.size();

    for (std::size_t i = 0; i < n; ++i) {
      const auto& edge = edges[i];

      if (query.index == edge.from_vertex) {
        idx = i; break;
      }
    }

    vs.insert(edges[(idx + n - 1) % n].from_vertex);
    vs.insert(edges[(idx + 1) % n].from_vertex);
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

  bool intersect_planes(
    const Plane_3& plane1,
    const Plane_3& plane2,
    Line_2& line_2) {

    auto result = CGAL::intersection(plane1, plane2);
    Line_3 line_3; bool found = false;
    if (result) {
      if (const Line_3* l = boost::get<Line_3>(&*result)) {
        found = true; line_3 = *l;
      }
    }
    if (!found) return false;

    const auto p1 = line_3.point(0);
    const auto p2 = line_3.point(1);
    const auto q1 = Point_2(p1.x(), p1.y());
    const auto q2 = Point_2(p2.x(), p2.y());

    line_2 = Line_2(q1, q2);
    return true;
  }

  /*
  void regularize_complex() {
    default_vertex_states();

    using Image_face_regularizer = internal::Image_face_complex_regularizer<
      Traits, Vertex, Edge, Halfedge, Face, Edge_type>;
    Image_face_regularizer image_face_regularizer(
      m_boundary,
      m_vertices, m_edges, m_halfedges,
      m_image, m_plane_map,
      m_min_length_2, m_angle_bound_2, m_ordinate_bound_2);

    std::vector<Edges> face_edges(m_faces.size());
    for (std::size_t i = 0; i < m_faces.size(); ++i) {
      const auto& face = m_faces[i];
      auto& edges = face_edges[i];
      create_face_edges(face, edges);
    }
    update_edge_neighbors(face_edges);

    for (std::size_t i = 0; i < m_faces.size(); ++i) {
      auto& face = m_faces[i];
      if (face.skip) continue;
      image_face_regularizer.set_face_edges(face_edges[i]);
      image_face_regularizer.compute_multiple_directions(face);
      image_face_regularizer.regularize_face(face);
    }

    // make_skip();

    default_vertex_states();
    for (auto& face : m_faces)
      update_face(face);

    save_faces_polylines("regularized");
    save_faces_ply("regularized");
  } */

  void default_vertex_states() {
    for (auto& vertex : m_vertices) {
      vertex.used  = false;
      vertex.state = false;
    }
  }

  void update_edge_neighbors(
    std::vector<Edges>& face_edges) {

    for (std::size_t i = 0; i < face_edges.size(); ++i) {
      auto& edges = face_edges[i];

      for (auto& edge : edges) {
        const std::size_t v1 = edge.from_vertex;
        const std::size_t v2 = edge.to_vertex;

        edge.faces.first  = i;
        edge.faces.second = find_neighbor_face(i, v1, v2, face_edges);
        if (edge.faces.second == std::size_t(-1))
          edge.type = Edge_type::BOUNDARY;
        else
          edge.type = Edge_type::INTERNAL;
      }
    }
  }

  std::size_t find_neighbor_face(
    const std::size_t skip,
    const std::size_t v1, const std::size_t v2,
    const std::vector<Edges>& face_edges) {

    for (std::size_t i = 0; i < face_edges.size(); ++i) {
      if (i == skip) continue;

      auto& edges = face_edges[i];
      for (auto& edge : edges) {
        const std::size_t w1 = edge.from_vertex;
        const std::size_t w2 = edge.to_vertex;

        if (v1 == w1 && v2 == w2)
          return i;
        if (v1 == w2 && v2 == w1)
          return i;
      }
    }
    return std::size_t(-1);
  }

  void make_skip() {
    for (auto& vertex : m_vertices)
      if (!vertex.skip)
        vertex.skip = vertex.state;
  }

  void initialize_vertices() {

    m_vertices.clear();
    for (const auto& ridge : m_ridges) {
      for (const auto& contour : ridge.contours) {
        if (contour.is_degenerated) continue;

        if (contour.is_closed)
          insert_vertices_closed(contour);
        else
          insert_vertices_open(contour);
      }
    }
    /* insert_vertices_boundary(); */
    clean_vertices();
    initialize_vertex_neighbors();
  }

  void clean_vertices() {

    std::size_t count = 0;
    std::vector<Vertex> clean;

    m_vertex_map.clear();
    for (std::size_t i = 0; i < m_vertices.size(); ++i) {
      auto& vertex = m_vertices[i];

      auto it = m_vertex_map.find(vertex.point);
      if (it == m_vertex_map.end()) {
        vertex.index = count; ++count;
        clean.push_back(vertex);
        m_vertex_map[vertex.point] = vertex.index;
      } else {
        auto& prev = clean[it->second];
        prev.type = vertex.type;
        prev.bd_idx = vertex.bd_idx;
        for (const std::size_t label : vertex.labels)
          prev.labels.insert(label);
      }
    }
    m_vertices = clean;
  }

  void insert_vertices_closed(
    const Contour& contour) {

    Vertex vertex;
    const auto& items = contour.points;
    const std::size_t m = items.size() - 1;

    for (std::size_t i = 0; i < m; ++i) {
      const auto& item = items[i];
      const auto& point = item.point();

      vertex.clear();
      vertex.point = point;
      vertex.labels = item.labels;
      vertex.type = item.end_type;
      m_vertices.push_back(vertex);
    }
  }

  void insert_vertices_open(
    const Contour& contour) {

    Vertex vertex;
    const auto& items = contour.points;
    const std::size_t m = items.size();
    CGAL_assertion(m >= 2);

    insert_corner_vertex(items[0]);
    for (std::size_t i = 1; i < m - 1; ++i) {
      const auto& item = items[i];
      const auto& point = item.point();

      vertex.clear();
      vertex.point = point;
      vertex.labels = item.labels;
      vertex.type = item.end_type;
      m_vertices.push_back(vertex);
    }
    insert_corner_vertex(items[m - 1]);
  }

  void insert_corner_vertex(const My_point& query) {

    for (auto& vertex : m_vertices) {
      if (internal::are_equal_points_2(vertex.point, query.point())) {
        for (const std::size_t label : query.labels)
          vertex.labels.insert(label);
        return;
      }
    }

    const bool is_boundary = ( query.end_type == Point_type::BOUNDARY );

    Vertex vertex;
    vertex.point = query.point();
    vertex.labels = query.labels;
    vertex.type = is_boundary ? Point_type::OUTER_BOUNDARY : Point_type::CORNER;
    vertex.bd_idx = is_boundary ? query.bd_idx : std::size_t(-1);
    m_vertices.push_back(vertex);
  }

  void insert_vertices_boundary() {

    Vertex vertex;
    for (std::size_t i = 0; i < m_boundary.size(); ++i) {
      const auto& segment = m_boundary[i];
      const auto& point = segment.source();

      for (auto& vertex : m_vertices) {
        if (internal::are_equal_points_2(vertex.point, point)) {

          vertex.type = Point_type::OUTER_CORNER;
          vertex.bd_idx = i;
          continue;
        }
      }

      vertex.clear();
      vertex.point = point;
      vertex.type = Point_type::OUTER_CORNER;
      vertex.bd_idx = i;
      m_vertices.push_back(vertex);
    }
  }

  void initialize_vertex_neighbors() {

    for (const auto& ridge : m_ridges) {
      for (const auto& contour : ridge.contours) {
        if (contour.is_degenerated) continue;

        if (contour.is_closed)
          insert_vertex_neighbors_closed(contour);
        else
          insert_vertex_neighbors_open(contour);
      }
    }
    /* insert_vertex_neighbors_boundary(); */
  }

  void insert_vertex_neighbors_closed(
    const Contour& contour) {

    if (contour.points.size() < 4) {
      std::cout.precision(30);
      std::cout << "Error: closed degenerated case!" << std::endl;
      for (const auto& item : contour.points)
        std::cout << item.point_ << std::endl;
      return;
    }

    const auto& items = contour.points;
    const std::size_t m = items.size() - 1;

    const auto& pr0 = items[m - 1].point();
    const auto& cr0 = items[0].point();
    const auto& nx0 = items[1].point();
    auto& ns0 = m_vertices[m_vertex_map.at(cr0)].neighbors;
    ns0.push_back(m_vertex_map.at(pr0));
    ns0.push_back(m_vertex_map.at(nx0));

    for (std::size_t i = 1; i < m - 1; ++i) {
      const std::size_t im = i - 1;
      const std::size_t ip = i + 1;

      const auto& prev = items[im].point();
      const auto& curr = items[i].point();
      const auto& next = items[ip].point();
      auto& neighbors = m_vertices[m_vertex_map.at(curr)].neighbors;
      neighbors.push_back(m_vertex_map.at(prev));
      neighbors.push_back(m_vertex_map.at(next));
    }

    const auto& pr1 = items[m - 2].point();
    const auto& cr1 = items[m - 1].point();
    const auto& nx1 = items[0].point();
    auto& ns1 = m_vertices[m_vertex_map.at(cr1)].neighbors;
    ns1.push_back(m_vertex_map.at(pr1));
    ns1.push_back(m_vertex_map.at(nx1));
  }

  void insert_vertex_neighbors_open(
    const Contour& contour) {

    if (contour.points.size() < 2) {
      std::cout.precision(30);
      std::cout << "Error: open degenerated case!" << std::endl;
      for (const auto& item : contour.points)
        std::cout << item.point_ << std::endl;
      return;
    }

    const auto& items = contour.points;
    const std::size_t m = items.size();

    insert_corner_neighbors(items[0], items[1]);
    for (std::size_t i = 1; i < m - 1; ++i) {
      const std::size_t im = i - 1;
      const std::size_t ip = i + 1;

      const auto& prev = items[im].point();
      const auto& curr = items[i].point();
      const auto& next = items[ip].point();
      auto& neighbors = m_vertices[m_vertex_map.at(curr)].neighbors;
      neighbors.push_back(m_vertex_map.at(prev));
      neighbors.push_back(m_vertex_map.at(next));
    }
    insert_corner_neighbors(items[m - 1], items[m - 2]);
  }

  void insert_corner_neighbors(
    const My_point& query, const My_point& other) {

    const auto& curr = query.point();
    auto& neighbors = m_vertices[m_vertex_map.at(curr)].neighbors;
    neighbors.push_back(m_vertex_map.at(other.point()));
  }

  void insert_vertex_neighbors_boundary() {

    std::vector<Indices> boundary_map;
    create_boundary_map(boundary_map);

    const std::size_t m = m_boundary.size();
    for (std::size_t i = 0; i < m; ++i) {

      const auto& curr = m_boundary[i].source();
      const std::size_t curr_idx = m_vertex_map.at(curr);
      auto& neighbors = m_vertices[curr_idx].neighbors;

      const std::size_t j = find_target_index(i, curr);
      CGAL_assertion(j != std::size_t(-1));

      add_one_boundary_point(curr, j, boundary_map, neighbors);
      add_one_boundary_point(curr, i, boundary_map, neighbors);
    }

    for (const auto& vertex : m_vertices) {
      if (vertex.type == Point_type::OUTER_BOUNDARY) {
        const auto& curr = vertex.point;
        const std::size_t curr_idx = m_vertex_map.at(curr);
        auto& neighbors = m_vertices[curr_idx].neighbors;
        add_two_boundary_points(
          curr_idx, vertex.bd_idx, boundary_map, neighbors);
      }
    }
  }

  void create_boundary_map(
    std::vector<Indices>& boundary_map) {

    const std::size_t m = m_boundary.size();

    boundary_map.clear();
    boundary_map.resize(m);

    for (std::size_t i = 0; i < m; ++i) {
      const auto& curr = m_boundary[i].source();
      const std::size_t curr_idx = m_vertex_map.at(curr);

      const auto& next = m_boundary[i].target();
      const std::size_t next_idx = m_vertex_map.at(next);

      boundary_map[i].push_back(curr_idx);
      boundary_map[i].push_back(next_idx);
    }

    for (const auto& vertex : m_vertices) {
      if (vertex.type == Point_type::OUTER_BOUNDARY) {
        const auto& curr = vertex.point;
        const std::size_t curr_idx = m_vertex_map.at(curr);
        boundary_map[vertex.bd_idx].push_back(curr_idx);
      }
    }
  }

  std::size_t find_target_index(
    const std::size_t skip, const Point_2& source) {

    const std::size_t m = m_boundary.size();
    for (std::size_t i = 0; i < m; ++i) {
      if (i == skip) continue;

      const auto& target = m_boundary[i].target();
      if (internal::are_equal_points_2(source, target))
        return i;
    }
    return std::size_t(-1);
  }

  void add_one_boundary_point(
    const Point_2& query,
    const std::size_t bd_idx,
    const std::vector<Indices>& boundary_map,
    Indices& neighbors) {

    auto ns = boundary_map[bd_idx];
    std::sort(ns.begin(), ns.end(),
    [this, &query](const std::size_t i, const std::size_t j) -> bool {
        const FT length_1 = internal::distance(query, m_vertices[i].point);
        const FT length_2 = internal::distance(query, m_vertices[j].point);
        return length_1 < length_2;
    });

    CGAL_assertion(ns.size() >= 2);
    neighbors.push_back(ns[1]);
  }

  void add_two_boundary_points(
    const std::size_t query_index,
    const std::size_t bd_idx,
    const std::vector<Indices>& boundary_map,
    Indices& neighbors) {

    const auto& ref = m_boundary[bd_idx].source();
    auto ns = boundary_map[bd_idx];

    std::sort(ns.begin(), ns.end(),
    [this, &ref](const std::size_t i, const std::size_t j) -> bool {
        const FT length_1 = internal::distance(ref, m_vertices[i].point);
        const FT length_2 = internal::distance(ref, m_vertices[j].point);
        return length_1 < length_2;
    });

    CGAL_assertion(ns.size() >= 3);
    for (std::size_t i = 1; i < ns.size() - 1; ++i) {
      const std::size_t im = i - 1;
      const std::size_t ip = i + 1;

      if (ns[i] == query_index) {
        neighbors.push_back(ns[im]);
        neighbors.push_back(ns[ip]);
        return;
      }
    }
  }

  void initialize_edges() {

    create_edges_and_halfedges();
    add_boundary_labels();
    set_edge_types();

    /* sort_vertex_halfedges(); */
    /* compute_next_halfedges(); */
  }

  void create_edges_and_halfedges() {

    m_edges.clear();
    m_halfedges.clear();
    m_halfedge_map.clear();

    std::size_t he_index = 0;
    std::size_t ed_index = 0;

    for (std::size_t i = 0; i < m_vertices.size(); ++i) {
      auto& vertexi = m_vertices[i];

      for (const std::size_t j : vertexi.neighbors) {
        auto& vertexj = m_vertices[j];

        CGAL_assertion(i != j);
        const auto to = std::make_pair(i, j);
        const auto op = std::make_pair(j, i);

        auto it = m_halfedge_map.find(op);
        if (it != m_halfedge_map.end()) {
          auto& other = it->second;

          Halfedge he;
          he.index = he_index; ++he_index;
          he.edg_idx = other.edg_idx;
          he.opposite = other.index;
          other.opposite = he.index;
          he.from_vertex = i;
          he.to_vertex = j;
          m_halfedge_map[to] = he;
          m_halfedges[other.index].opposite = he.index;
          m_halfedges.push_back(he);
          update_edge_labels(
            vertexi, vertexj, m_edges[he.edg_idx]);

          vertexi.hedges.push_back(he.index);

        } else {

          Edge edge;
          edge.index = ed_index; ++ed_index;
          update_edge_labels(
            vertexi, vertexj, edge);
          edge.from_vertex = i;
          edge.to_vertex = j;
          edge.segment = Segment_2(vertexi.point, vertexj.point);
          m_edges.push_back(edge);

          Halfedge he;
          he.index = he_index; ++he_index;
          he.edg_idx = edge.index;
          he.from_vertex = i;
          he.to_vertex = j;
          m_halfedge_map[to] = he;
          m_halfedges.push_back(he);

          vertexi.hedges.push_back(he.index);
        }
      }
    }
    m_halfedge_map.clear();
  }

  void update_edge_labels(
    const Vertex& vertexi, const Vertex& vertexj,
    Edge& edge) {

    if (
      edge.labels.first  != std::size_t(-1) &&
      edge.labels.second != std::size_t(-1) )
    return;

    if (
      vertexi.type != Point_type::OUTER_BOUNDARY &&
      vertexi.type != Point_type::OUTER_CORNER) {

      if (vertexi.labels.size() <= 2) {
        auto it = vertexi.labels.begin();

        if (vertexi.labels.size() >= 1) {
          edge.labels.first = *it;
        }
        if (vertexi.labels.size() == 2) {
          ++it; edge.labels.second = *it;
        }

      } else if (vertexj.labels.size() <= 2) {
        auto it = vertexj.labels.begin();

        if (vertexj.labels.size() >= 1) {
          edge.labels.first = *it;
        }
        if (vertexj.labels.size() == 2) {
          ++it; edge.labels.second = *it;
        }
      } else {

        bool foundi = false;
        for (const std::size_t label : vertexi.labels) {
          if (label == std::size_t(-1)) {
            foundi = true; break;
          }
        }
        bool foundj = false;
        for (const std::size_t label : vertexj.labels) {
          if (label == std::size_t(-1)) {
            foundj = true; break;
          }
        }
        const bool is_boundary = foundi && foundj;

        if (is_boundary)
          add_boundary_label(vertexi, vertexj, edge);
        else
          add_internal_label(vertexi, vertexj, edge);

        /*
        std::cout << "s:" <<
        vertexi.labels.size() << " " << vertexj.labels.size() << std::endl;
        if (
          edge.labels.first == std::size_t(-1) ||
          edge.labels.second == std::size_t(-1)) {

          std::cout.precision(30);
          std::cout << vertexi.point << std::endl;
          std::cout << vertexj.point << std::endl;
        }
        std::cout <<
        edge.labels.first << " " << edge.labels.second << std::endl; */
      }
    }
  }

  void add_internal_label(
    const Vertex& vertexi, const Vertex& vertexj,
    Edge& edge) {

    const auto& s = vertexi.point;
    const auto& t = vertexj.point;
    const auto  m = internal::middle_point_2(s, t);
    set_internal_labels(m, edge);
  }

  void set_internal_labels(
    const Point_2& query, Edge& edge) {

    Indices neighbors;
    m_knq(query, neighbors);

    auto& labels = edge.labels;
    const auto& pixels = m_image.pixels;
    for (const std::size_t n : neighbors) {
      if (
        pixels[n].label != std::size_t(-1) &&
        labels.first == std::size_t(-1)) {

        labels.first = pixels[n].label;
        continue;
      }
      if (
        labels.first != std::size_t(-1) &&
        pixels[n].label != labels.first) {

        labels.second = pixels[n].label;
        break;
      }
    }
  }

  void add_boundary_label(
    const Vertex& vertexi, const Vertex& vertexj,
    Edge& edge) {

    const auto& s = vertexi.point;
    const auto& t = vertexj.point;
    const auto  m = internal::middle_point_2(s, t);
    set_boundary_labels(m, edge);
  }

  void set_boundary_labels(
    const Point_2& query, Edge& edge) {

    Indices neighbors;
    m_knq(query, neighbors);

    auto& labels = edge.labels;
    const auto& pixels = m_image.pixels;
    for (const std::size_t n : neighbors) {
      if (pixels[n].label != std::size_t(-1)) {
        labels.first = pixels[n].label;
        return;
      }
    }
  }

  void add_boundary_labels() {

    for (auto& edge : m_edges) {
      auto& labels = edge.labels;
      if (
        labels.first  == std::size_t(-1) &&
        labels.second == std::size_t(-1)) {
        add_boundary_label(edge);
      }
    }
  }

  void add_boundary_label(Edge& edge) {

    const auto& s = m_vertices[edge.from_vertex].point;
    const auto& t = m_vertices[edge.to_vertex].point;
    const auto  m = internal::middle_point_2(s, t);
    set_labels(m, edge);
  }

  void set_labels(
    const Point_2& query, Edge& edge) {

    Indices neighbors;
    m_knq(query, neighbors);

    auto& labels = edge.labels;
    const auto& pixels = m_image.pixels;
    for (const std::size_t n : neighbors) {
      if (pixels[n].label != std::size_t(-1)) {
        labels.first = pixels[n].label;
        return;
      }
    }
  }

  void set_edge_types() {

    for (auto& edge : m_edges) {
      const auto& labels = edge.labels;
      edge.type = Edge_type::DEFAULT;

      const std::size_t l1 = labels.first;
      const std::size_t l2 = labels.second;

      if (l1 == std::size_t(-1) || l2 == std::size_t(-1)) {
        edge.type = Edge_type::BOUNDARY; continue;
      }

      const auto& from = m_vertices[edge.from_vertex];
      const auto& to   = m_vertices[edge.to_vertex];

      if (from.type != Point_type::LINEAR || to.type != Point_type::LINEAR) {
        edge.type = Edge_type::INTERNAL; continue;
      }
    }
  }

  void sort_vertex_halfedges() {

    for (auto& vertex : m_vertices) {
      auto& hedges = vertex.hedges;

      std::sort(hedges.begin(), hedges.end(),
      [this](const std::size_t i, const std::size_t j) -> bool {

        const auto& hedgei = m_halfedges[i];
        const auto& hedgej = m_halfedges[j];

        const std::size_t idx_toi = hedgei.to_vertex;
        const std::size_t idx_toj = hedgej.to_vertex;

        const auto typei = m_vertices[idx_toi].type;
        const auto typej = m_vertices[idx_toj].type;

        const std::size_t pi = get_priority(typei);
        const std::size_t pj = get_priority(typej);

        return pi > pj;
      });
    }
  }

  std::size_t get_priority(const Point_type type) {

    switch (type) {
      case Point_type::FREE:
        return 4;
      case Point_type::LINEAR:
        return 4;
      case Point_type::CORNER:
        return 3;
      case Point_type::OUTER_BOUNDARY:
        return 2;
      case Point_type::OUTER_CORNER:
        return 1;
      default:
        return 0;
    }
    return 0;
  }

  void compute_next_halfedges() {

    for (const auto& he : m_halfedges) {
      if (he.next != std::size_t(-1)) continue;
      if (m_halfedges[he.opposite].next != std::size_t(-1)) continue;

      const auto& labels = m_edges[he.edg_idx].labels;

      const std::size_t l1 = labels.first;
      const std::size_t l2 = labels.second;
      if (l1 == std::size_t(-1) || l2 == std::size_t(-1))
        continue;

      traverse(he.index, l1);
      traverse(he.opposite, l2);
    }
  }

  void traverse(
    const std::size_t idx,
    const std::size_t ref_label,
    Indices& hedges) {

    hedges.clear();
    const auto& he = m_halfedges[idx];
    if (he.next != std::size_t(-1)) return;
    if (ref_label == std::size_t(-1)) return;
    if (m_halfedges[he.opposite].next != std::size_t(-1)) return;

    std::vector<Segment_2> segments;

    std::size_t count = 0;
    const std::size_t start = he.index;
    std::size_t curr = start;
    do {

      hedges.push_back(curr);
      auto& other = m_halfedges[curr];
      const std::size_t to_idx = other.to_vertex;
      const auto& to = m_vertices[to_idx];
      find_next(start, ref_label, to, other);

      const auto& s = m_vertices[other.from_vertex].point;
      const auto& t = m_vertices[other.to_vertex].point;

      segments.push_back(Segment_2(s, t));
      curr = other.next;

      /*
      if (m_halfedges[curr].next != std::size_t(-1)) {
        other.next = std::size_t(-1); return;
      }
      if (m_halfedges[m_halfedges[curr].opposite].next != std::size_t(-1)) {
        other.next = std::size_t(-1); return;
      } */

      if (count >= 10000) {

        std::cout.precision(30);
        std::cout << "Error: traverse() max count reached!" << std::endl;
        std::cout << "Ref label: " << ref_label << std::endl;

        Saver saver;
        saver.save_polylines(
        segments, "/Users/monet/Documents/gf/lod/logs/buildings/tmp/debug-count");
        saver.save_polylines(
        m_boundary, "/Users/monet/Documents/gf/lod/logs/buildings/tmp/debug-boundary-count");
        exit(EXIT_FAILURE);
      }

      if (curr == std::size_t(-1)) {
        std::cout << "Error: traverse() failed!" << std::endl;
        std::cout << "Ref label: " << ref_label << std::endl;

        Saver saver;
        saver.save_polylines(
        segments, "/Users/monet/Documents/gf/lod/logs/buildings/tmp/debug-fail");
        saver.save_polylines(
        m_boundary, "/Users/monet/Documents/gf/lod/logs/buildings/tmp/debug-boundary-fail");
        exit(EXIT_FAILURE);
      }
      ++count;

    } while (curr != start);
  }

  void find_next(
    const std::size_t start,
    const std::size_t ref_label,
    const Vertex& to, Halfedge& he) {

    for (const std::size_t other_idx : to.hedges) {
      const auto& other = m_halfedges[other_idx];
      if (other.edg_idx == he.edg_idx)
        continue;

      if (other.index == start) {
        he.next = other.index; return;
      }

      if (
        other.next != std::size_t(-1) ||
        m_halfedges[other.opposite].next != std::size_t(-1)) {
        return;
      }

      const auto& edge = m_edges[other.edg_idx];
      const auto& labels = edge.labels;
      const std::size_t l1 = labels.first;
      const std::size_t l2 = labels.second;

      if (l1 == ref_label && l2 != ref_label) {
        he.next = other.index;
        /* traverse(other.opposite, l2); */
        return;
      }
      if (l2 == ref_label && l1 != ref_label) {
        he.next = other.index;
        /* traverse(other.opposite, l1); */
        return;
      }
      if (l1 != ref_label && l2 != ref_label) {
        continue;
      }

      std::cout << "Error: find_next() failed!" << std::endl;
      exit(EXIT_FAILURE);
    }
  }

  void initialize_faces() {

    m_faces.clear();
    std::vector<Indices> regions;

    /*
    regions.clear();
    add_face_regions(regions);
    create_faces(regions); */

    std::set<std::size_t> labels;
    for (const auto& edge : m_edges) {

      if (edge.labels.first != std::size_t(-1))
        labels.insert(edge.labels.first);
      if (edge.labels.second != std::size_t(-1))
        labels.insert(edge.labels.second);
    }

    for (const std::size_t ref_label : labels) {
      /* compute_next_he(ref_label); */
      compute_next_vt(ref_label, regions);
      /*
      regions.clear();
      add_face_regions(regions);
      add_faces(ref_label, regions); */
    }

    if (m_faces.size() == 0) {
      const std::size_t ref_label = *labels.begin();
      compute_next_he(ref_label);

      /*
      regions.clear();
      add_face_regions(regions);
      add_faces(ref_label, regions); */
    }

    clear_next_halfedges();
    remove_duplicated_faces();
    for (std::size_t i = 0; i < m_faces.size(); ++i) {
      m_faces[i].index = i;
      initialize_face(m_faces[i]);
    }
    clean_faces();

    sort_faces();
    create_face_neighbors();
    set_face_types();

    for (auto& face : m_faces) {
      for (auto fh = face.tri.delaunay.finite_faces_begin();
      fh != face.tri.delaunay.finite_faces_end(); ++fh)
        fh->info().label = face.label;
    }
  }

  void clean_faces() {

    std::vector<Face> clean;
    for (auto& face : m_faces) {
      if (!is_contained(face))
        clean.push_back(face);
    }
    m_faces.clear();
    m_faces.reserve(clean.size());
    for (std::size_t i = 0; i < clean.size(); ++i) {
      auto& face = clean[i];
      face.index = i;
      m_faces.push_back(face);
    }
  }

  bool is_contained(const Face& query) {

    FT max_area = -FT(1);
    Point_2 best;

    const auto& delaunay = query.tri.delaunay;
    for (auto fh = delaunay.finite_faces_begin();
    fh != delaunay.finite_faces_end(); ++fh) {
      if (!fh->info().tagged) continue;

      const auto& p0 = fh->vertex(0)->point();
      const auto& p1 = fh->vertex(1)->point();
      const auto& p2 = fh->vertex(2)->point();

      const Point_2 b = CGAL::barycenter(
        p0, FT(1), p1, FT(1), p2, FT(1));

      const Triangle_2 triangle = Triangle_2(p0, p1, p2);
      const FT area = CGAL::abs(triangle.area());
      if (area > max_area) {
        max_area = area; best = b;
      }
    }

    if (max_area != -FT(1)) {
      for (const auto& face : m_faces) {
        if (query.index == face.index) continue;
        const auto& del = face.tri.delaunay;
        const auto fh = del.locate(best);
        if (fh->info().tagged) {

          if (face.area > query.area)
            return true;
          else return false;
        }
      }
    } else return true;
    return false;
  }

  void set_face_types() {

    for (auto& face : m_faces) {

      std::size_t count = 0;
      for (const std::size_t he_idx : face.hedges) {
        const auto& vertex = m_vertices[m_halfedges[he_idx].from_vertex];
        for (const std::size_t label : vertex.labels) {
          if (label == std::size_t(-1)) {
            ++count; break;
          }
        }
      }

      if (count > 1)
        face.type = Face_type::BOUNDARY;
      else
        face.type = Face_type::INTERNAL;
    }
  }

  void compute_next_he(
    const std::size_t ref_label) {

    clear_next_halfedges();
    Face face;
    for (const auto& he : m_halfedges) {
      if (he.next != std::size_t(-1)) continue;
      if (m_halfedges[he.opposite].next != std::size_t(-1)) continue;

      const auto& edge = m_edges[he.edg_idx];
      const auto& labels = edge.labels;

      const std::size_t l1 = labels.first;
      const std::size_t l2 = labels.second;

      /*
      if (l1 == std::size_t(-1) || l2 == std::size_t(-1))
        continue; */

      if (l1 == ref_label) {
        clear_next_halfedges();
        traverse(he.index, l1, face.hedges);
        face.label = ref_label;
        m_faces.push_back(face);
        clear_next_halfedges();
        return;
      }

      if (l2 == ref_label) {
        clear_next_halfedges();
        traverse(he.index, l2, face.hedges);
        face.label = ref_label;
        m_faces.push_back(face);
        clear_next_halfedges();
        return;
      }
    }
  }

  void compute_next_vt(
    const std::size_t ref_label,
    std::vector<Indices>& regions) {

    regions.clear();
    clear_next_halfedges();

    Face face;
    for (const auto& vt : m_vertices) {
      if (!(
        vt.type == Point_type::CORNER ||
        vt.type == Point_type::OUTER_BOUNDARY ||
        vt.type == Point_type::OUTER_CORNER) ) continue;

      for (const std::size_t he_idx : vt.hedges) {
        const auto& he = m_halfedges[he_idx];

        if (he.next != std::size_t(-1)) continue;
        if (m_halfedges[he.opposite].next != std::size_t(-1)) continue;

        const auto& edge = m_edges[he.edg_idx];
        const auto& labels = edge.labels;

        const std::size_t l1 = labels.first;
        const std::size_t l2 = labels.second;

        /*
        if (l1 == std::size_t(-1) || l2 == std::size_t(-1))
          continue; */

        if (l1 == ref_label || l2 == ref_label) {
          clear_next_halfedges();
          traverse(he_idx, ref_label, face.hedges);
          if (face.hedges.size() != 0) {
            face.label = ref_label;
            m_faces.push_back(face);
          }
          /*
          add_face_regions(regions);
          add_faces(ref_label, regions);
          regions.clear(); */
          clear_next_halfedges();
        }
      }
    }
  }

  void clear_next_halfedges() {
    for (auto& he : m_halfedges)
      he.next = std::size_t(-1);
  }

  void add_face_regions(
    std::vector<Indices>& regions) {

    DS_neighbor_query neighbor_query(
      m_vertices, m_edges, m_halfedges);
    DS_region_type region_type(
      m_vertices, m_edges, m_halfedges);
    Region_growing region_growing(
      m_halfedges, neighbor_query, region_type);

    region_growing.detect(std::back_inserter(regions));
    /* std::cout << "num face regions: " << regions.size() << std::endl; */
  }

  void add_faces(
    const std::size_t ref_label,
    const std::vector<Indices>& regions) {

    if (regions.size() == 0)
      return;

    Face face;
    for (const auto& region : regions) {
      face.hedges = region;
      face.label = ref_label;
      m_faces.push_back(face);
    }
  }

  /*
  void create_faces(
    const std::vector<Indices>& regions) {

    const std::size_t numr = regions.size();
    if (numr == 0) return;

    m_faces.clear();
    m_faces.reserve(numr);

    Face face;
    for (std::size_t i = 0; i < numr - 1; ++i) {
      if (are_equal_regions(regions[i], regions[i + 1]))
        continue;

      face.hedges.clear();
      face.hedges = regions[i];
      m_faces.push_back(face);
    }

    face.hedges.clear();
    face.hedges = regions[numr - 1];
    m_faces.push_back(face);
  }
  */

  bool are_equal_regions(
    const Indices& r1, const Indices& r2) {

    if (r1.size() != r2.size())
      return false;

    const std::size_t numr = r1.size();
    for (std::size_t i = 0; i < numr; ++i) {
      const std::size_t j = numr - i - 1;

      const std::size_t ei = m_halfedges[r1[i]].edg_idx;
      const std::size_t ej = m_halfedges[r2[j]].edg_idx;
      if (ei != ej) return false;
    }
    return true;
  }

  void remove_duplicated_faces() {

    std::vector<Face> clean;
    std::vector<bool> states(m_faces.size(), false);

    Indices dindices;
    for (std::size_t i = 0; i < m_faces.size(); ++i) {
      if (states[i]) continue;

      const auto& f1 = m_faces[i];
      find_duplicates(i, f1, states, dindices);
      for (const std::size_t didx : dindices)
        states[didx] = true;
    }

    for (std::size_t i = 0; i < m_faces.size(); ++i)
      if (!states[i]) clean.push_back(m_faces[i]);
    m_faces = clean;
  }

  void find_duplicates(
    const std::size_t skip,
    const Face& f1,
    const std::vector<bool>& states,
    Indices& dindices) {

    dindices.clear();
    for (std::size_t i = 0; i < m_faces.size(); ++i) {
      if (i == skip || states[i]) continue;

      const auto& f2 = m_faces[i];
      if (are_equal_faces(f1, f2))
        dindices.push_back(i);
    }
  }

  bool are_equal_faces(
    const Face& f1, const Face& f2) {

    const auto& hedges1 = f1.hedges;
    const auto& hedges2 = f2.hedges;

    if (f1.label != f2.label)
      return false;

    if (hedges1.size() != hedges2.size())
      return false;

    std::size_t count = 0;
    for (const std::size_t idx1 : hedges1) {
      const std::size_t edg_idx1 = m_halfedges[idx1].edg_idx;
      for (const std::size_t idx2 : hedges2) {
        const std::size_t edg_idx2 = m_halfedges[idx2].edg_idx;
        if (edg_idx1 == edg_idx2) {
          ++count; break;
        }
      }
    }

    /*
    std::cout <<
    f1.label << " " << f2.label << " " <<
    hedges1.size() << " " << hedges2.size() <<
    " " << count << std::endl; */

    return ( hedges1.size() == count );
  }

  void initialize_face(Face& face) {

    face.type = Face_type::DEFAULT;
    create_face_probs(face);
    create_face_triangulation(face);
    create_face_visibility(face);
    /* create_face_label(face); */
    compute_face_area(face);
  }

  void create_face_neighbors() {

    for (auto& vertex : m_vertices)
      vertex.tmp.clear();

    for (auto& face : m_faces) {
      face.used = false;
      face.tmp.clear();
    }

    for (std::size_t i = 0; i < m_faces.size(); ++i) {
      auto& facei = m_faces[i];
      for (std::size_t j = 0; j < m_faces.size(); ++j) {
        if (j == i) continue;
        auto& facej = m_faces[j];
        if (!facej.used)
          add_face_neighbors(facei, facej);
      }
    }

    for (auto& vertex : m_vertices) {
      vertex.faces.clear();
      vertex.faces.reserve(vertex.tmp.size());
      for (const std::size_t fidx : vertex.tmp)
        vertex.faces.push_back(fidx);
      vertex.tmp.clear();
    }

    for (auto& face : m_faces) {
      face.used = false;
      face.neighbors.clear();
      face.neighbors.reserve(face.tmp.size());

      for (const std::size_t fidx : face.tmp)
        face.neighbors.push_back(fidx);
      face.tmp.clear();
    }
  }

  void add_face_neighbors(
    Face& f1, Face& f2) {

    for (const std::size_t h1 : f1.hedges) {
      auto& v1 = m_vertices[m_halfedges[h1].from_vertex];
      for(const std::size_t h2 : f2.hedges) {
        auto& v2 = m_vertices[m_halfedges[h2].from_vertex];

        if (v1.index == v2.index) {
          v1.tmp.insert(f1.index);
          v1.tmp.insert(f2.index);
          f1.tmp.insert(f2.index);
          f2.tmp.insert(f1.index);
        }
      }
    }
  }

  bool update_face(Face& face) {

    if (face.skip) return false;
    const bool success = create_face_triangulation(face);
    if (!success) {
      std::cout << "Warning: Empty face, update_face(), Image_data_structure!" << std::endl;
      return false;
    }
    create_face_visibility(face);
    update_face_label(face);
    compute_face_area(face);
    if (face.area < internal::tolerance<FT>())
      face.skip = true;
    return true;
  }

  void create_face_probs(Face& face) {

    face.probs.clear();
    for (const std::size_t he_idx : face.hedges) {
      const auto& he = m_halfedges[he_idx];
      const auto& edge = m_edges[he.edg_idx];
      const auto& labels = edge.labels;

      if (labels.first != std::size_t(-1))
        face.probs.insert(labels.first);
      if (labels.second != std::size_t(-1))
        face.probs.insert(labels.second);
    }
  }

  bool create_face_triangulation(Face& face) {

    auto& tri = face.tri.delaunay;
    tri.clear();

    std::vector<Edge> edges;
    create_face_edges(face, edges);

    if (edges.size() == 0)
      return false;

    for (const auto& edge : edges) {
      const auto& s = edge.segment.source();
      const auto& t = edge.segment.target();

      const auto vh1 = tri.insert(s);
      const auto vh2 = tri.insert(t);
      vh1->info().object_index = edge.from_vertex;
      vh2->info().object_index = edge.to_vertex;

      if (vh1 != vh2)
        tri.insert_constraint(vh1, vh2);
    }

    if (tri.number_of_faces() == 0)
      return false;
    return true;
  }

  void create_face_visibility(
    Face& face) {

    std::vector<Point_2> bbox;
    create_face_bbox(face, bbox);

    auto& tri = face.tri.delaunay;
    for (auto fh = tri.finite_faces_begin();
    fh != tri.finite_faces_end(); ++fh) {
      fh->info().object_index = face.index;

      const auto& p0 = fh->vertex(0)->point();
      const auto& p1 = fh->vertex(1)->point();
      const auto& p2 = fh->vertex(2)->point();

      const FT x = (p0.x() + p1.x() + p2.x()) / FT(3);
      const FT y = (p0.y() + p1.y() + p2.y()) / FT(3);
      const Point_2 p = Point_2(x, y);

      FT in = FT(1); FT out = FT(1);
      for (std::size_t i = 0; i < bbox.size(); ++i) {
        const auto& q = bbox[i];
        if (tri.oriented_side(fh, p) == CGAL::ON_NEGATIVE_SIDE)
          continue;

        LF_circulator circ = tri.line_walk(p, q, fh);
        const LF_circulator end = circ;
        if (circ.is_empty()) continue;

        std::size_t inter = 0;
        do {

          LF_circulator f1 = circ; ++circ;
          LF_circulator f2 = circ;

          const bool success = are_neighbors(f1, f2);
          if (!success) break;

          const std::size_t idx = f1->index(f2);
          const auto edge = std::make_pair(f1, idx);
          if (tri.is_constrained(edge)) ++inter;
          if (tri.is_infinite(f2)) break;

        } while (circ != end);

        if (inter % 2 == 0) out += FT(1);
        else in += FT(1);
      }

      const FT sum = in + out;
      in /= sum; out /= sum;

      if (in > FT(1) / FT(2)) {
        fh->info().interior = true;
        fh->info().tagged   = true;
      } else {
        fh->info().interior = false;
        fh->info().tagged   = false;
      }
    }
  }

  void create_face_bbox(
    const Face& face,
    std::vector<Point_2>& bbox) {

    bbox.clear();
    const auto& tri = face.tri.delaunay;

    std::vector<Point_2> points;
    for (auto fh = tri.finite_faces_begin();
    fh != tri.finite_faces_end(); ++fh) {
      const auto& p0 = fh->vertex(0)->point();
      const auto& p1 = fh->vertex(1)->point();
      const auto& p2 = fh->vertex(2)->point();
      points.push_back(p0);
      points.push_back(p1);
      points.push_back(p2);
    }

    CGAL::Identity_property_map<Point_2> pmap;
    internal::bounding_box_2(points, pmap, bbox);
  }

  bool are_neighbors(
    LF_circulator f1, LF_circulator f2) const {

    for (std::size_t i = 0; i < 3; ++i) {
      const std::size_t ip = (i + 1) % 3;

      const auto p1 = f1->vertex(i);
      const auto p2 = f1->vertex(ip);

      for (std::size_t j = 0; j < 3; ++j) {
        const std::size_t jp = (j + 1) % 3;

        const auto q1 = f2->vertex(j);
        const auto q2 = f2->vertex(jp);

        if (
          ( p1 == q1 && p2 == q2) ||
          ( p1 == q2 && p2 == q1) ) {

          return true;
        }
      }
    }
    return false;
  }

  void create_face_label(Face& face) {
    clear_face_probabilities(face);
    compute_face_probabilities(face);
    initialize_face_labels(face);
    set_face_label(face);
  }

  void clear_face_probabilities(Face& face) {

    const std::size_t num_labels = m_image.num_labels;
    auto& tri = face.tri.delaunay;
    for (auto fh = tri.finite_faces_begin();
    fh != tri.finite_faces_end(); ++fh) {
      if (!fh->info().interior) continue;

      fh->info().probabilities.clear();
      fh->info().probabilities.resize(num_labels, FT(0));
    }
  }

  void compute_face_probabilities(Face& face) {

    auto& tri = face.tri.delaunay;
    for (const auto& pixel : m_image.pixels) {
      if (pixel.label == std::size_t(-1))
        continue;

      const auto& p = pixel.point;
      auto fh = tri.locate(p);
      if (tri.is_infinite(fh) || !fh->info().interior)
        continue;
      fh->info().probabilities[pixel.label] += FT(1);
    }
  }

  void initialize_face_labels(Face& face) {

    const std::size_t num_labels = m_image.num_labels;
    auto& tri = face.tri.delaunay;
    for (auto fh = tri.finite_faces_begin();
    fh != tri.finite_faces_end(); ++fh) {
      if (!fh->info().interior) continue;

      FT max_prob = FT(-1);
      std::size_t best_label = std::size_t(-1);
      for (std::size_t i = 0; i < num_labels; ++i) {
        const FT prob = fh->info().probabilities[i];
        if (prob > max_prob) {
          max_prob = prob; best_label = i;
        }
      }

      if (max_prob != FT(0))
        fh->info().label = best_label;
      else
        fh->info().label = std::size_t(-1);
    }
  }

  void set_face_label(Face& face) {

    const std::size_t num_labels = m_image.num_labels;
    std::vector<std::size_t> probs(num_labels, FT(0));

    auto& tri = face.tri.delaunay;
    for (auto fh = tri.finite_faces_begin();
    fh != tri.finite_faces_end(); ++fh) {
      if (!fh->info().interior) continue;
      if (fh->info().label != std::size_t(-1))
        probs[fh->info().label] += FT(1);
    }

    FT max_prob = FT(-1);
    std::size_t best_label = std::size_t(-1);
    for (std::size_t i = 0; i < num_labels; ++i) {
      if (probs[i] > max_prob) {
        max_prob = probs[i]; best_label = i;
      }
    }

    for (auto fh = tri.finite_faces_begin();
    fh != tri.finite_faces_end(); ++fh) {
      if (!fh->info().interior) {
        fh->info().label = std::size_t(-1); continue;
      }
      fh->info().label = best_label;
    }
    face.label = best_label;
  }

  void update_face_label(Face& face) {

    auto& tri = face.tri.delaunay;
    for (auto fh = tri.finite_faces_begin();
    fh != tri.finite_faces_end(); ++fh) {
      if (!fh->info().interior) {
        fh->info().label = std::size_t(-1); continue;
      }
      fh->info().label = face.label;
    }
  }

  void compute_face_area(Face& face) {

    face.area = FT(0);
    const auto& tri = face.tri.delaunay;
    for (auto fh = tri.finite_faces_begin();
    fh != tri.finite_faces_end(); ++fh) {
      if (!fh->info().interior) continue;

      const auto& p0 = fh->vertex(0)->point();
      const auto& p1 = fh->vertex(1)->point();
      const auto& p2 = fh->vertex(2)->point();

      const Triangle_2 triangle = Triangle_2(p0, p1, p2);
      face.area += triangle.area();
    }
  }

  void sort_faces() {

    std::sort(m_faces.begin(), m_faces.end(),
    [](const Face& f1, const Face& f2) -> bool {
      return f1.area > f2.area;
    });
    for (std::size_t i = 0; i < m_faces.size(); ++i)
      m_faces[i].index = i;
  }

  void mark_bad_faces() {

    /* mark_bad_faces_tri_based(); */
    mark_bad_faces_area_based();

    for (const auto& face : m_faces) {
      if (face.skip) {

        bool corner_found = false, bound_found = false;
        FT x = FT(0), y = FT(0), count = FT(0);
        for (const std::size_t he_idx : face.hedges) {
          const auto& he = m_halfedges[he_idx];
          auto& vertex = m_vertices[he.from_vertex];
          if (
            vertex.type == Point_type::LINEAR ||
            vertex.type == Point_type::FREE) {

            vertex.skip = true;
          }

          if (!corner_found && vertex.type == Point_type::OUTER_BOUNDARY) {
            x += vertex.point.x();
            y += vertex.point.y();
            count += FT(1);
            bound_found = true;
          }
          if (!corner_found && vertex.type == Point_type::OUTER_CORNER) {
            x = vertex.point.x();
            y = vertex.point.y();
            corner_found = true;
          }
        }

        if (!corner_found && bound_found) {
          x /= count; y /= count;
        }

        if (corner_found || bound_found) {
          for (const std::size_t he_idx : face.hedges) {
            const auto& he = m_halfedges[he_idx];
            auto& vertex = m_vertices[he.from_vertex];

            if (
              vertex.type == Point_type::OUTER_BOUNDARY ||
              vertex.type == Point_type::OUTER_CORNER) {

              vertex.point = Point_2(x, y);
            }
          }
        }
      }
    }
    for (auto& face : m_faces) {
      const bool success = update_face(face);
      if (!success) face.skip = true;
    }
  }

  void mark_bad_faces_tri_based() {
    for (auto& face : m_faces)
      if (face.tri.delaunay.number_of_faces() == 1)
        face.skip = true;
  }

  void mark_bad_faces_area_based() {

    /*
    FT avg_area = FT(0), count = FT(0);
    for (const auto& face : m_faces) {
      if (face.skip) continue;
      avg_area += face.area;
      count += FT(1);
    }
    avg_area /= count; */

    for (auto& face : m_faces) {
      if (face.skip) continue;
      if (face.area < FT(1) / FT(2))
        face.skip = true;
    }
  }

  void create_face_edges(
    const Face& face,
    std::vector<Edge>& edges) {

    edges.clear();
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

        edges.push_back(edge);
        continue;
      }

      if (!s.skip && t.skip) {
        bool last = false;
        i = get_next(face, i, last);

        const std::size_t next_idx = face.hedges[i];
        const auto& next = m_halfedges[next_idx];
        const std::size_t end = next.to_vertex;
        const auto& other = m_vertices[end];

        Edge edge;
        edge.from_vertex = from;
        edge.to_vertex = end;
        edge.segment = Segment_2(s.point, other.point);
        edge.type = Edge_type::INTERNAL;

        edges.push_back(edge);
        if (last) break;
        continue;
      }
    }
  }

  std::size_t get_next(
    const Face& face,
    const std::size_t start,
    bool& last) {

    for (std::size_t i = start; i < face.hedges.size(); ++i) {
      const std::size_t he_idx = face.hedges[i];
      const auto& he = m_halfedges[he_idx];
      const std::size_t to = he.to_vertex;
      if (m_vertices[to].skip) continue;
      return i;
    }

    last = true;
    for (std::size_t i = 0; i < start; ++i) {
      const std::size_t he_idx = face.hedges[i];
      const auto& he = m_halfedges[he_idx];
      const std::size_t to = he.to_vertex;
      if (m_vertices[to].skip) continue;
      return i;
    }

    return std::size_t(-1);
  }

  Point_3 get_position_on_plane_3(
    const std::size_t vidx,
    const Point_2& query,
    const Plane_3& ref) {

    const auto refp = internal::position_on_plane_3(query, ref);
    if (vidx == std::size_t(-1))
      return refp;

    const auto& vertex = m_vertices[vidx];
    const auto& faces = vertex.faces;

    FT z = FT(0); FT count = FT(0);
    for (const std::size_t fidx : faces) {
      if (fidx == std::size_t(-1)) continue;
      const auto& face = m_faces[fidx];
      if (face.skip) continue;
      const std::size_t label = face.label;
      if (label == std::size_t(-1)) continue;

      const auto& plane = m_plane_map.at(label);
      const FT val = internal::position_on_plane_3(query, plane).z();
      if (CGAL::abs(val - refp.z()) < m_max_height_difference / FT(2)) {
        z += val; count += FT(1);
      }
    }

    if (count == FT(0))
      return refp;
    z /= count;
    return Point_3(query.x(), query.y(), z);
  }

  void snap_to_boundary() {

    std::vector<Segment_2> oriented;
    orient_boundary(oriented);
    m_boundary = oriented;
    for (std::size_t i = 0; i < m_boundary.size(); ++i) {
      const auto& segment = m_boundary[i];
      const auto& query = segment.source();
      set_closest_corner(i, query);
    }
    set_outer_boundary_points();
    project_boundaries();
  }

  void orient_boundary(
    std::vector<Segment_2>& contour) {

    Saver saver;
    contour.clear();
    bool finished = false;
    std::size_t count_out = 0;
    std::vector<bool> states(m_boundary.size(), false);
    m_ends.clear();
    std::size_t init = 0;

    do {
      bool stop_out = false;
      auto curr  = find_longest_segment(m_boundary, states, stop_out);
      if (stop_out) {
        saver.save_polylines(
          contour, "/Users/monet/Documents/gf/lod/logs/buildings/tmp/stop-out");
        std::cout << "Warning: stop out reached!" << std::endl;
        break;
      }
      auto start = curr;

      finished = true;
      bool completed = false;
      std::size_t count_in = 0;

      do {
        contour.push_back(curr);

        bool stop_in = false;
        curr = find_next_curr(curr, states, stop_in);
        if (stop_in) {
          saver.save_polylines(
          contour, "/Users/monet/Documents/gf/lod/logs/buildings/tmp/stop-in");
          std::cout << "Warning: stop in reached!" << std::endl;
          break;
        }

        ++count_in;
        if (curr.target() == start.source()) {
          contour.push_back(curr);
          for (std::size_t i = init; i < contour.size(); ++i)
            m_ends[i] = std::make_pair(init, contour.size());
          init = contour.size();
          completed = true;
        }
      } while (!completed && count_in != 10000);

      if (count_in == 10000) {
        std::cout << "Error: max count in reached, orient_boundary()!" << std::endl;

        saver.save_polylines(
        contour, "/Users/monet/Documents/gf/lod/logs/buildings/tmp/debug-count");
        exit(EXIT_FAILURE);
      }
      ++count_out;

      for (const bool state : states)
        if (!state) finished = false;
    } while (!finished && count_out != 100);

    if (count_out == 100) {
      std::cout << "Error: max count out reached, orient_boundary()!" << std::endl;
      exit(EXIT_FAILURE);
    }
  }

  Segment_2 find_longest_segment(
    const std::vector<Segment_2>& segments,
    std::vector<bool>& states,
    bool& stop) {

    FT max_length = -FT(1);
    std::size_t seg_idx = std::size_t(-1);

    for (std::size_t i = 0; i < segments.size(); ++i) {
      if (states[i]) continue;
      const auto& segment = segments[i];
      const FT length = segment.squared_length();
      if (length > max_length) {

        max_length = length;
        seg_idx = i;
      }
    }

    if (seg_idx != std::size_t(-1)) {
      states[seg_idx] = true;
      return segments[seg_idx];
    }
    stop = true;
    return Segment_2();
  }

  Segment_2 find_next_curr(
    const Segment_2& curr,
    std::vector<bool>& states,
    bool& stop) {

    for (std::size_t i = 0; i < m_boundary.size(); ++i) {
      if (states[i]) continue;

      const auto& segment = m_boundary[i];
      if (
        internal::are_equal_points_2(curr.source(), segment.source()) &&
        internal::are_equal_points_2(curr.target(), segment.target()))
        continue;

      if (
        internal::are_equal_points_2(curr.source(), segment.target()) &&
        internal::are_equal_points_2(curr.target(), segment.source()))
        continue;

      if (internal::are_equal_points_2(curr.target(), segment.source())) {
        states[i] = true;
        return segment;
      }

      if (internal::are_equal_points_2(curr.target(), segment.target())) {
        states[i] = true;
        return Segment_2(segment.target(), segment.source());
      }
    }
    stop = true;
    return Segment_2();
  }

  void set_closest_corner(
    const std::size_t bd_idx,
    const Point_2& query) {

    FT min_dist = internal::max_value<FT>();
    std::size_t closest = std::size_t(-1);
    for (std::size_t i = 0; i < m_vertices.size(); ++i) {
      const auto& vertex = m_vertices[i];
      if (vertex.used) continue;

      bool is_boundary = false;
      for (const std::size_t label : vertex.labels) {
        if (label == std::size_t(-1)) {
          is_boundary = true; break;
        }
      }
      if (!is_boundary) continue;

      const FT dist = internal::distance(query, vertex.point);
      if (dist < min_dist) {
        min_dist = dist;
        closest = i;
      }
    }
    m_vertices[closest].point  = query;
    m_vertices[closest].used   = true;
    m_vertices[closest].type   = Point_type::OUTER_CORNER;
    m_vertices[closest].bd_idx = bd_idx;
  }

  void set_outer_boundary_points() {

    for (auto& vertex : m_vertices) {
      if (vertex.used) continue;
      if (vertex.type == Point_type::OUTER_CORNER)
        continue;

      bool is_boundary = false;
      for (const std::size_t label : vertex.labels) {
        if (label == std::size_t(-1)) {
          is_boundary = true; break;
        }
      }
      if (!is_boundary) continue;

      vertex.type = Point_type::BOUNDARY;
      if (vertex.labels.size() > 2)
        vertex.type = Point_type::OUTER_BOUNDARY;
    }
  }

  void project_boundaries() {
    for (std::size_t i = 0; i < m_vertices.size(); ++i) {
      auto& vertex = m_vertices[i];
      if (vertex.type == Point_type::BOUNDARY)
        vertex.skip = true;
    }
    project_outer_boundaries();
  }

  void project_outer_boundaries() {

    Indices neighbors;
    for (std::size_t i = 0; i < m_vertices.size(); ++i) {
      auto& vertex = m_vertices[i];
      if (vertex.type == Point_type::OUTER_BOUNDARY) {

        neighbors.clear();
        for (const std::size_t neighbor : vertex.neighbors)
          if (is_boundary_neighbor(neighbor))
            neighbors.push_back(neighbor);
        CGAL_assertion(neighbors.size() >= 2);

        const std::size_t c1 = find_corner(i, neighbors[0]);
        const std::size_t c2 = find_corner(i, neighbors[1]);
        CGAL_assertion(c1 != c2);

        if (c1 == std::size_t(-1) || c2 == std::size_t(-1))
          continue;

        const std::size_t bd_idx1 = m_vertices[c1].bd_idx;
        const std::size_t bd_idx2 = m_vertices[c2].bd_idx;

        const std::size_t b1 = m_ends.at(bd_idx1).first;
        const std::size_t b2 = m_ends.at(bd_idx1).second;

        if (bd_idx1 == b1 && bd_idx2 == b2 - 1) {
          project_onto_boundary(c2, vertex); continue;
        }
        if (bd_idx1 == b2 - 1 && bd_idx2 == b1) {
          project_onto_boundary(c1, vertex); continue;
        }

        if (bd_idx1 < bd_idx2)
          project_onto_boundary(c1, vertex);
        else
          project_onto_boundary(c2, vertex);
      }
    }
  }

  bool is_boundary_neighbor(const std::size_t idx) {
    bool is_boundary = false;
    for (const std::size_t label : m_vertices[idx].labels)
      if (label == std::size_t(-1))
        return true;
    return false;
  }

  std::size_t find_corner(
    const std::size_t start,
    const std::size_t idx) {

    std::size_t prev = start;
    std::size_t curr = idx;
    bool is_corner = false;

    std::size_t count = 0;
    do {
      const auto& vertex = m_vertices[curr];
      if (vertex.type == Point_type::OUTER_CORNER)
        return vertex.index;

      bool found = false;
      for (const std::size_t neighbor : vertex.neighbors) {
        if (neighbor == prev) continue;
        if (!is_boundary_neighbor(neighbor)) continue;
        if (m_vertices[neighbor].type == Point_type::OUTER_CORNER)
          return neighbor;

        found = true; prev = curr;
        curr = neighbor; break;
      }

      if (!found) {
        std::cerr <<
          "Error: neighbor not found, find_corner()!"
          << std::endl;
        exit(EXIT_FAILURE);
      }
      ++count;
    } while (!is_corner && count <= 10000);

    if (count >= 10000) {
      std::cerr <<
      "Error: neighbor not found, max count reached find_corner()!"
      << std::endl;
    }
    return std::size_t(-1);
  }

  void project_onto_boundary(
    const std::size_t corner_idx,
    Vertex& vertex) {

    const std::size_t bd_idx = m_vertices[corner_idx].bd_idx;
    const auto& segment = m_boundary[bd_idx];
    const Line_2 line = Line_2(segment.source(), segment.target());
    auto proj = line.projection(vertex.point);
    vertex.point = proj;
    vertex.bd_idx = bd_idx;
  }
};

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_INTERNAL_IMAGE_DATA_STRUCTURE_H
