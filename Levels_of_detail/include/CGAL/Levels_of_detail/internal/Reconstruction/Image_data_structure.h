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
#include <CGAL/property_map.h>

// Internal includes.
#include <CGAL/Levels_of_detail/internal/utils.h>
#include <CGAL/Levels_of_detail/internal/struct.h>

// Other includes.
#include <CGAL/Levels_of_detail/internal/Reconstruction/Image.h>

// Testing.
#include "../../../../../test/Levels_of_detail/include/Saver.h"

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

template<typename GeomTraits>
struct Image_data_structure {

public:
  using Traits = GeomTraits;

  using FT = typename Traits::FT;
  using Point_2 = typename Traits::Point_2;
  using Segment_2 = typename Traits::Segment_2;
  
  using Image = internal::Image<Traits>;
  
  using My_point = typename Image::My_point;
  using Contour = typename Image::Contour;
  using Point_type = typename Image::Point_type;

  using Indices = std::vector<std::size_t>;
  using Size_pair = std::pair<std::size_t, std::size_t>;

  enum class Face_type {
    DEFAULT = 0,
    CLOSED = 1,
    BASIC = 2,
    INTERNAL = 3
  };

  struct Vertex {
    Point_2 point;
    std::set<std::size_t> labels;
    std::size_t index = std::size_t(-1);
    std::size_t bd_idx = std::size_t(-1);
    Point_type type = Point_type::DEFAULT;
    Indices edges;
    Indices neighbors;

    void clear() {
      labels.clear();
      point = Point_2();
      index = std::size_t(-1);
      bd_idx = std::size_t(-1);
      type = Point_type::DEFAULT;
      edges.clear();
      neighbors.clear();
    }
  };

  struct Edge {
    Size_pair labels = std::make_pair(std::size_t(-1), std::size_t(-1));
    Segment_2 segment;
    std::size_t index = std::size_t(-1);
  };

  struct Halfedge {
    std::size_t index = std::size_t(-1);
    std::size_t edg_idx = std::size_t(-1);
    std::size_t next = std::size_t(-1);
    std::size_t opposite = std::size_t(-1);
  };

  struct Face {
    Indices edges;
    Indices neighbors;
    std::size_t label = std::size_t(-1);
    std::size_t type = Face_type::DEFAULT;
  };

  Image_data_structure(
    const std::vector<Segment_2>& boundary,
    const std::vector<Image>& ridges) :
  m_boundary(boundary),
  m_ridges(ridges),
  m_pi(static_cast<FT>(CGAL_PI)) 
  { }

  void build() {

    initialize_vertices();
    initialize_vertex_neighbors();

    /*
    for (const auto& vertex : m_vertices) {
      if (
        vertex.type != Point_type::FREE && 
        vertex.type != Point_type::LINEAR) {

        std::cout << 
        int(vertex.type) << " : " << 
        vertex.bd_idx << " , " <<
        vertex.labels.size() << " , " << 
        vertex.neighbors.size() << std::endl;
      }
    } */

    initialize_edges();

    std::cout << "num edges: " << m_edges.size() << std::endl;
    std::cout << "num halfedges: " << m_halfedges.size() << std::endl;
    for (const auto& he : m_halfedges) {
      std::cout << he.second.edg_idx << " : " << he.second.opposite << std::endl;
    }
  }

  void clear() {
    m_vertices.clear();
    m_edges.clear();
    m_faces.clear();
    m_halfedges.clear();
  }

private:
  const std::vector<Segment_2>& m_boundary;
  const std::vector<Image>& m_ridges;
  const FT m_pi;

  std::vector<Vertex> m_vertices;
  std::vector<Edge>   m_edges;
  std::vector<Face>   m_faces;

  std::map<Point_2, std::size_t> m_vertex_map;
  std::map<Size_pair, Halfedge> m_halfedges;

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
    insert_vertices_boundary();
    
    m_vertex_map.clear();
    for (std::size_t i = 0; i < m_vertices.size(); ++i) {
      m_vertices[i].index = i;
      m_vertex_map[m_vertices[i].point] = m_vertices[i].index;
    }
  }

  void insert_vertices_closed(
    const Contour& contour) {

    Vertex vertex;
    const auto& items = contour.points;
    const std::size_t m = items.size() - 1;

    const std::size_t start = m_vertices.size();
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
    const std::size_t start = m_vertices.size();
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
    insert_vertex_neighbors_boundary();
  }

  void insert_vertex_neighbors_closed(
    const Contour& contour) {

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

    std::map<std::size_t, Indices> boundary_map;
    create_boundary_map(boundary_map);

    const std::size_t m = m_boundary.size();
    for (std::size_t i = 0; i < m; ++i) {
      const std::size_t im = (i + m - 1) % m;

      const auto& curr = m_boundary[i].source();
      const std::size_t curr_idx = m_vertex_map.at(curr);
      auto& neighbors = m_vertices[curr_idx].neighbors;
      add_one_boundary_point(curr, im, boundary_map, neighbors);
      add_one_boundary_point(curr, i , boundary_map, neighbors);
    }

    for (const auto& vertex : m_vertices) {
      if (vertex.type == Point_type::OUTER_BOUNDARY) {
        const auto& curr = vertex.point;
        const std::size_t curr_idx = m_vertex_map.at(curr);
        auto& neighbors = m_vertices[curr_idx].neighbors;
        add_two_boundary_points(curr, vertex.bd_idx, boundary_map, neighbors);
      }
    }
  }

  void create_boundary_map(
    std::map<std::size_t, Indices>& boundary_map) {
    
    boundary_map.clear();
    const std::size_t m = m_boundary.size();

    for (std::size_t i = 0; i < m; ++i) {
      const std::size_t ip = (i + 1) % m;

      const auto& curr = m_boundary[i].source();
      const std::size_t curr_idx = m_vertex_map.at(curr);

      const auto& next = m_boundary[ip].source();
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

  void add_one_boundary_point(
    const Point_2& query,
    const std::size_t bd_idx, 
    const std::map<std::size_t, Indices>& boundary_map,
    Indices& neighbors) {

    auto ns = boundary_map.at(bd_idx);
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
    const Point_2& query,
    const std::size_t bd_idx, 
    const std::map<std::size_t, Indices>& boundary_map,
    Indices& neighbors) {

    auto ns = boundary_map.at(bd_idx);
    std::sort(ns.begin(), ns.end(), 
    [this, &query](const std::size_t i, const std::size_t j) -> bool { 
        const FT length_1 = internal::distance(query, m_vertices[i].point);
        const FT length_2 = internal::distance(query, m_vertices[j].point);
        return length_1 < length_2;
    });

    CGAL_assertion(ns.size() >= 3);
    neighbors.push_back(ns[1]);
    neighbors.push_back(ns[2]);
  }

  void initialize_edges() {

    std::size_t he_index = 0;
    std::size_t ed_index = 0;

    for (std::size_t i = 0; i < m_vertices.size(); ++i) {
      auto& vertexi = m_vertices[i];

      for (const std::size_t j : vertexi.neighbors) {
        auto& vertexj = m_vertices[j];

        const auto to = std::make_pair(i, j);
        const auto op = std::make_pair(j, i);

        auto it = m_halfedges.find(op);
        if (it != m_halfedges.end()) {
          auto& other = it->second;

          Halfedge he;
          he.index = he_index; ++he_index;
          he.edg_idx = other.edg_idx;
          he.opposite = other.index;
          other.opposite = he.index;
          m_halfedges[to] = he;

        } else {
          
          Edge edge;
          edge.index = ed_index; ++ed_index;

          if (vertexi.labels.size() <= 2) {
            auto it = vertexi.labels.begin();

            if (vertexi.labels.size() >= 1) {
              edge.labels.first = *it;
            }
            if (vertexi.labels.size() >= 2) {
              ++it; edge.labels.first = *it;
            }
            
          } else if (vertexj.labels.size()) {
            auto it = vertexj.labels.begin();

            if (vertexj.labels.size() >= 1) {
              edge.labels.first = *it;
            }
            if (vertexj.labels.size() >= 2) {
              ++it; edge.labels.first = *it;
            }
          }

          edge.segment = Segment_2(vertexi.point, vertexj.point);
          m_edges.push_back(edge);

          Halfedge he;
          he.index = he_index; ++he_index;
          he.edg_idx = edge.index;
          m_halfedges[to] = he;
        }
      }
    }
  }
};

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_INTERNAL_IMAGE_DATA_STRUCTURE_H
