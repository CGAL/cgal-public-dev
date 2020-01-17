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

#ifndef CGAL_LEVELS_OF_DETAIL_INTERNAL_IMAGE_TREE_H
#define CGAL_LEVELS_OF_DETAIL_INTERNAL_IMAGE_TREE_H

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
class Image_tree {

public:
  using Traits = GeomTraits;

  using FT = typename Traits::FT;
  using Point_2 = typename Traits::Point_2;
  using Point_3 = typename Traits::Point_3;
  using Vector_2 = typename Traits::Vector_2;
  using Segment_2 = typename Traits::Segment_2;
  using Segment_3 = typename Traits::Segment_3;
  
  using Image = internal::Image<Traits>;
  using Point_type = typename Image::Point_type;

  using FT_pair = std::pair<FT, FT>;
  using Indices = std::vector<std::size_t>;
  using Size_pair = std::pair<std::size_t, std::size_t>;

  using Saver = Saver<Traits>;
  using Inserter = Polygon_inserter<Traits>;
  using Indexer = internal::Indexer<Point_3>;
  using Color = CGAL::Color;

  using Edges = std::vector<Edge>;

  enum class Node_type {
    DEFAULT = 0,
    ROOT = 1,
    CHILD = 2,
    LEAF = 3
  };

  struct Node {
    std::size_t index = std::size_t(-1);
    Node_type type = Node_type::DEFAULT;
    std::size_t face_index = std::size_t(-1);
    std::size_t label = std::size_t(-1);
    Indices children;

    void add_children(
      const std::vector<Node>& nodes,
      Indices& faces) const {

      for (const std::size_t child : children) {
        const auto& node = nodes[child];

        faces.push_back(node.face_index);
        node.add_children(nodes, faces);
      }
    }

    void clear() {
      children.clear();
    }
  };

  struct Tree {
    std::vector<Node> nodes;
    std::vector<Indices> levels;

    void traverse_children(
      const std::size_t node_idx,
      Indices& faces) {
      
      faces.clear();
      const auto& node = nodes[node_idx];
      node.add_children(nodes, faces);
    }

    void clear() {
      nodes.clear();
      levels.clear();
    }
  };

  Image_tree(
    const std::vector<Segment_2>& boundary,
    std::vector<Vertex>& vertices,
    std::vector<Edge>& edges,
    std::vector<Halfedge>& halfedges, 
    std::vector<Face>& faces) :
  m_boundary(boundary),
  m_vertices(vertices),
  m_edges(edges),
  m_halfedges(halfedges),
  m_faces(faces),
  m_pi(static_cast<FT>(CGAL_PI)),
  m_bound_min(FT(15)),
  m_bound_max(FT(75)) { 

    m_directions.clear();
    const std::size_t seg_idx = find_longest_segment(m_boundary);
    m_directions.push_back(m_boundary[seg_idx]);

    std::sort(m_directions.begin(), m_directions.end(), 
    [](const Segment_2& a, const Segment_2& b) -> bool { 
      return a.squared_length() > b.squared_length();
    });
  }

  void build() {  
    build_tree();
  }

  void cut(std::size_t level) {

    if (level < 0) level = 0;
    if (level >= m_tree.levels.size())
      level = m_tree.levels.size() - 1;
    cut_along_tree(level);
  }

  void check_vertex_information() {

    std::cout.precision(30);
    std::cout << "num vertices: " << m_vertices.size() << std::endl;
    for (const auto& vertex : m_vertices) {
      if (
        vertex.type != Point_type::FREE && 
        vertex.type != Point_type::LINEAR) {
        
        std::cout << 
        int(vertex.type) << " : " << 
        vertex.bd_idx << " , " <<
        vertex.labels.size() << " , " << 
        vertex.hedges.size() << " , " <<
        vertex.neighbors.size() << " , " <<
        vertex.faces.size() << std::endl;
      }
    }
  }

  void check_edge_information() {
    
    std::cout << "num edges: " << m_edges.size() << std::endl;
    std::cout << "num halfedges: " << m_halfedges.size() << std::endl;
    for (const auto& edge : m_edges) {
      std::cout << 
        int(edge.type) << " : " << 
        edge.labels.first << " " << 
        edge.labels.second << std::endl;
    }
  }

  void check_halfedge_information() {
    
    std::cout << "num edges: " << m_edges.size() << std::endl;
    std::cout << "num halfedges: " << m_halfedges.size() << std::endl;
    for (const auto& he : m_halfedges) {
      std::cout <<
      he.from_vertex << " " << he.to_vertex << " : " <<  
      int(m_vertices[he.from_vertex].type) << " " <<
      int(m_vertices[he.to_vertex].type) << " : " <<
      he.index << " " << 
      he.edg_idx << " " << 
      he.next << " " << 
      he.opposite << std::endl;
    }
  }

  void check_face_information() {
    
    std::cout << "num faces: " << m_faces.size() << std::endl;
    for (const auto& face : m_faces) {
      std::cout <<
      face.index << " : " << 
      int(face.type) << " , " <<
      face.label << " , " <<
      face.probs.size() << " , " <<
      face.hedges.size() << " , " <<
      face.neighbors.size() << " , " <<
      face.area << " , " <<
      face.skip << std::endl;
    }
  }

  void check_tree_information(
    const bool check_nodes = true,
    const bool check_levels = true) {

    if (check_nodes) {
      std::cout << "num nodes: " << m_tree.nodes.size() << std::endl;
      for (const auto& node : m_tree.nodes) {
        std::cout << 
        node.index << " : " << 
        int(node.type) << " , " << 
        node.face_index << " , " <<
        node.label << " , " <<
        node.children.size() << std::endl;
      }
    }

    if (check_levels) {
      std::cout << "num levels: " << m_tree.levels.size() << std::endl;
      for (std::size_t i = 0; i < m_tree.levels.size(); ++i)
        std::cout << i << " : " << m_tree.levels[i].size() << std::endl;
    }
  }

  std::size_t num_levels() {
    return m_tree.levels.size();
  }

private:
  const std::vector<Segment_2>& m_boundary;
  std::vector<Vertex>& m_vertices;
  std::vector<Edge>& m_edges;
  std::vector<Halfedge>& m_halfedges;
  std::vector<Face>& m_faces;
  
  const FT m_pi;
  const FT m_bound_min;
  const FT m_bound_max;
  
  Tree m_tree;
  std::vector<Segment_2> m_directions;

  std::size_t find_longest_segment(
    const std::vector<Segment_2>& segments) {

    std::size_t seg_idx = std::size_t(-1);
    FT max_length = -FT(1);
    for (std::size_t i = 0; i < segments.size(); ++i) {
        
      const FT length = segments[i].squared_length();
      if (length > max_length) {

        max_length = length;
        seg_idx = i;
      }
    }
    return seg_idx;
  }

  void build_tree() {
    m_tree.clear();
    create_tree_nodes();
    create_tree_levels();
  }

  void create_tree_nodes() {
    CGAL_assertion(m_faces.size() != 0);

    auto& nodes = m_tree.nodes;
    nodes.clear();
    nodes.reserve(m_faces.size());

    Node node;
    node.index = 0;

    if (m_faces.size() == 1) {
      node.face_index = 0;
      node.label = m_faces[node.face_index].label;
      m_faces[node.face_index].level = 0;
    } else {
      node.face_index = std::size_t(-1);
      node.label = std::size_t(-1);
    }
    
    node.type = Node_type::ROOT;
    nodes.push_back(node);

    if (m_faces.size() > 1) {
      for (std::size_t i = 0; i < m_faces.size(); ++i) {
        const auto& face = m_faces[i];

        node.index = i + 1;
        node.face_index = face.index;
        node.label = face.label;
        node.type = Node_type::CHILD;
        nodes.push_back(node);
      }
    }
  }

  void create_tree_levels() {

    m_tree.levels.clear();
    std::vector<Edges> face_edges(m_faces.size());
    for (std::size_t i = 0; i < m_faces.size(); ++i) {
      const auto& face = m_faces[i];
      auto& edges = face_edges[i];
      create_face_edges(face, edges);
    }
    update_edge_neighbors(face_edges);

    create_root_level();
    create_base_level(face_edges);
    create_mansard_level(face_edges);
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

  void create_root_level() {
    Indices root_indices(1, 0);
    m_tree.levels.push_back(root_indices);
  }

  void create_base_level(
    const std::vector<Edges>& face_edges) {

    auto& nodes  = m_tree.nodes;
    auto& levels = m_tree.levels;
    
    const FT avg_area = compute_average_face_area();
    const FT eps = avg_area / FT(2);

    Indices level;
    for (std::size_t i = 0; i < face_edges.size(); ++i) {
      const auto& edges = face_edges[i];
      auto& face = m_faces[i];
      if (face.area < eps) continue;

      for (const auto& edge : edges) {
        if (comply_with_outer_boundary(edge)) {
          level.push_back(i + 1); 
          face.level = 1; break;
        }
      }
    }

    levels.push_back(level);
    auto& root = nodes[0];
    root.children = level;
  }

  bool comply_with_outer_boundary(
    const Edge& edge) {

    const auto& segment = edge.segment;
    const FT slength = 
      internal::distance(segment.source(), segment.target());

    for (const auto& direction : m_directions) {
      const FT dlength = 
        internal::distance(direction.source(), direction.target());
      
      if (slength < dlength / FT(2))
        continue;

      const FT angle   = angle_degree_2(direction, segment);
      const FT angle_2 = get_angle_2(angle);

      if ( 
        (CGAL::abs(angle_2) <= m_bound_min) ||
        (CGAL::abs(angle_2) >= m_bound_max) )  {
        return true;
      }
    }
    return false;
  }

  FT compute_average_face_area() {

    CGAL_assertion(m_faces.size() != 0);
    FT avg_area = FT(0);
    for (const auto& face : m_faces)
      avg_area += face.area;
    avg_area /= static_cast<FT>(m_faces.size());
    return avg_area;
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

  void create_mansard_level(
    const std::vector<Edges>& face_edges) {

    auto& nodes  = m_tree.nodes;
    auto& levels = m_tree.levels;

    Indices level;
    for (std::size_t i = 0; i < face_edges.size(); ++i) {
      const auto& edges = face_edges[i];
      auto& face = m_faces[i];
      if (face.level != std::size_t(-1)) continue;
      
      const auto& neighbors = face.neighbors;
      const std::size_t best_face = get_best_face(edges, neighbors);
      if (best_face == std::size_t(-1)) continue;

      const std::size_t node_idx = best_face + 1;
      level.push_back(i + 1);
      nodes[node_idx].children.push_back(i + 1);
      face.level = 2;
    }
    levels.push_back(level);
  }

  std::size_t get_best_face(
    const Edges& edges,
    const Indices& neighbors) {

    std::map<std::size_t, FT> length;
    for (const std::size_t idx : neighbors) {
      if (m_faces[idx].level != 1) continue;
      length[idx] = FT(0);
    }

    for (const auto& edge : edges) {
      const std::size_t f2 = edge.faces.second;
      const auto& seg = edge.segment;

      if (length.find(f2) != length.end())
        length[f2] += internal::distance(seg.source(), seg.target());
    }

    FT max_dist = -FT(1);
    std::size_t best_idx = std::size_t(-1);

    for (const auto& pair : length) {
      if (pair.second != FT(0)) {
        if (max_dist < pair.second) {
          max_dist = pair.second;
          best_idx = pair.first;
        }
      }
    }
    return best_idx;
  }

  void cut_along_tree(const std::size_t lidx) {
    
    Indices faces;
    const auto& nodes = m_tree.nodes;
    const auto& level = m_tree.levels[lidx];

    if (nodes.size() == 1) {
      m_faces[0].label = nodes[0].label; return;
    }

    for (std::size_t i = 1; i < nodes.size(); ++i)
      m_faces[i - 1].label = nodes[i].label;

    for (const std::size_t nidx : level) {
      m_tree.traverse_children(nidx, faces);
      for (const std::size_t fidx : faces) {
        auto& face = m_faces[fidx];
        face.label = nodes[nidx].label;
      }
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

        edges.push_back(edge);
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
};

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_INTERNAL_IMAGE_TREE_H
