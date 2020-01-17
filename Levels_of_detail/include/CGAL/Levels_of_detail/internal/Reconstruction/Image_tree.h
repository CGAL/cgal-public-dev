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
    std::vector<Vertex>& vertices,
    std::vector<Edge>& edges,
    std::vector<Halfedge>& halfedges, 
    std::vector<Face>& faces) :
  m_vertices(vertices),
  m_edges(edges),
  m_halfedges(halfedges),
  m_faces(faces),
  m_pi(static_cast<FT>(CGAL_PI)) 
  { }

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
  std::vector<Vertex>& m_vertices;
  std::vector<Edge>& m_edges;
  std::vector<Halfedge>& m_halfedges;
  std::vector<Face>& m_faces;
  const FT m_pi;
  Tree m_tree;

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

    Node node; std::size_t count = 0;
    node.index = count; ++count;

    if (m_faces.size() == 1) {
      node.face_index = 0;
      node.label = m_faces[node.face_index].label;
    } else {
      node.face_index = std::size_t(-1);
      node.label = std::size_t(-1);
    }
    
    node.type = Node_type::ROOT;
    nodes.push_back(node);

    if (m_faces.size() > 1) {
      for (const auto& face : m_faces) {
        node.index = count; ++count;
        node.face_index = face.index;
        node.label = m_faces[node.face_index].label;
        node.type = Node_type::CHILD;
        nodes.push_back(node);
      }
    }
  }

  void create_tree_levels() {

    auto& nodes = m_tree.nodes;
    auto& levels = m_tree.levels;
    
    // Root.
    levels.clear();
    Indices root_indices(1, 0);
    levels.push_back(root_indices);
    
    auto& root = nodes[0];
    root.children.clear();

    // Update face boundaries.
    std::vector< std::vector<Edge> > face_edges(m_faces.size());
    for (std::size_t i = 0; i < m_faces.size(); ++i) {
      const auto& face = m_faces[i];
      auto& edges = face_edges[i];
      create_face_edges(face, edges);
    }

    // Create children.
    /*
    for (std::size_t i = 1; i < nodes.size(); ++i)
      root.children.push_back(i); */
  }

  void cut_along_tree(const std::size_t lidx) {
    
    Indices faces;
    const auto& nodes = m_tree.nodes;
    const auto& level = m_tree.levels[lidx];

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
