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

// Boost includes.
#include <boost/shared_ptr.hpp>

// CGAL includes.
#include <CGAL/enum.h>
#include <CGAL/property_map.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/create_straight_skeleton_2.h>

#define CGAL_DO_NOT_USE_BOYKOV_KOLMOGOROV_MAXFLOW_SOFTWARE
#include <CGAL/internal/Surface_mesh_segmentation/Alpha_expansion_graph_cut.h>

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
  using Plane_3 = typename Traits::Plane_3;
  
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
  using Alpha_expansion = CGAL::internal::Alpha_expansion_graph_cut_boost;

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

    void print_children(
      const std::vector<Node>& nodes) const {
      for (const std::size_t child : children)
        std::cout << child << " ";
      std::cout << std::endl;
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
    const std::vector<Segment_2>& directions,
    const std::map<std::size_t, Plane_3>& plane_map,
    const FT min_length_2,
    const FT max_height_difference,
    std::vector<Vertex>& vertices,
    std::vector<Edge>& edges,
    std::vector<Halfedge>& halfedges, 
    std::vector<Face>& faces) :
  m_boundary(boundary),
  m_directions(directions),
  m_plane_map(plane_map),
  m_min_length_2(min_length_2),
  m_max_height_difference(max_height_difference),
  m_vertices(vertices),
  m_edges(edges),
  m_halfedges(halfedges),
  m_faces(faces),
  m_pi(static_cast<FT>(CGAL_PI)),
  m_bound_min(FT(10)),
  m_bound_max(FT(80)),
  m_beta(FT(1) / FT(2)) { 

    /*
    m_directions.clear();
    const std::size_t seg_idx = find_longest_segment(m_boundary);
    m_directions.push_back(m_boundary[seg_idx]); */

    std::cout << "num directions: " << m_directions.size() << std::endl;
    std::sort(m_directions.begin(), m_directions.end(), 
    [](const Segment_2& a, const Segment_2& b) -> bool { 
      return a.squared_length() > b.squared_length();
    });
  }

  void build() {
    build_tree_with_graphcut();
  }

  void build_v1() {
    build_skeleton();
    build_tree_naive();
  }

  void apply_test() {

    for (auto& face : m_faces) {
      if (face.neighbors.size() == 1)
        face.label = m_faces[face.neighbors[0]].label;
      
      for (auto fh = face.tri.delaunay.finite_faces_begin(); 
      fh != face.tri.delaunay.finite_faces_end(); ++fh)
        fh->info().label = face.label;
    }

    Indices indices;
    indices.reserve(m_faces.size());
    for (std::size_t i = 0; i < m_faces.size(); ++i) 
      indices.push_back(i);
    
    std::sort(indices.begin(), indices.end(), 
    [&](const std::size_t i, const std::size_t j) -> bool { 
      return m_faces[i].neighbors.size() < m_faces[j].neighbors.size();
    });
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
  const std::map<std::size_t, Plane_3>& m_plane_map;
  
  const FT m_min_length_2;
  const FT m_max_height_difference;

  std::vector<Vertex>& m_vertices;
  std::vector<Edge>& m_edges;
  std::vector<Halfedge>& m_halfedges;
  std::vector<Face>& m_faces;
  
  const FT m_pi;
  const FT m_bound_min;
  const FT m_bound_max;
  const FT m_beta;
  
  Tree m_tree;
  std::vector<Segment_2> m_directions;

  std::map<std::size_t, std::size_t> m_dr_mapping;
  std::map<std::size_t, std::size_t> m_op_mapping;

  void build_skeleton() {
    
    using Polygon_2 = CGAL::Polygon_2<Traits>;

    Polygon_2 polygon;
    for (const auto& segment : m_boundary)
      polygon.push_back(segment.source());
    auto iss = CGAL::create_interior_straight_skeleton_2(polygon);

    std::vector<Segment_2> segments;
    for (auto he = iss->halfedges_begin(); 
    he != iss->halfedges_end(); ++he) {
      if (
        he->is_bisector() && 
        ( ( he->id() % 2 ) == 0 ) && 
        !he->has_infinite_time() && 
        !he->opposite()->has_infinite_time()) { 

        auto segment = Segment_2(
          he->vertex()->point(), he->opposite()->vertex()->point());
        if (has_boundary_vertex(segment)) continue;
        segments.push_back(segment);
      }
    }

    Saver saver;
    saver.save_polylines(
    segments, "/Users/monet/Documents/lod/logs/buildings/tmp/skeleton");
  }

  bool has_boundary_vertex(
    const Segment_2& query) {

    for (const auto& segment : m_boundary) {
      if (query.source() == segment.source()) return true;
      if (query.target() == segment.source()) return true;
      if (query.source() == segment.target()) return true;
      if (query.target() == segment.target()) return true;
    }
    return false;
  }

  void print_tree() {

    const auto& nodes  = m_tree.nodes;
    const auto& levels = m_tree.levels;

    std::cout << std::endl;
    for (std::size_t i = 0; i < levels.size(); ++i) {
      std::cout << "level " << i << " : " << std::endl;

      for (const std::size_t idx : levels[i]) {
        std::cout << "node " << idx << " has children" << " -> " << std::endl;
        nodes[idx].print_children(nodes);
      }
      std::cout << std::endl;
    }
  }

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

  void build_tree_with_graphcut() {
    m_tree.clear();
    create_tree_nodes();
    create_tree_levels_with_graphcut();
  }

  void build_tree_naive() {
    m_tree.clear();
    create_tree_nodes();
    create_tree_levels_naive();
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

  void create_tree_levels_naive() {

    m_tree.levels.clear();
    std::vector<Edges> face_edges(m_faces.size());
    for (std::size_t i = 0; i < m_faces.size(); ++i) {
      const auto& face = m_faces[i];
      auto& edges = face_edges[i];
      create_face_edges(face, edges);
    }
    update_edge_neighbors(face_edges);

    create_root_level();
    create_base_level_naive(face_edges);
    create_sublevels_naive(face_edges);
  }

  void create_tree_levels_with_graphcut() {
    
    m_tree.levels.clear();
    compute_face_weights();
    std::cout << "faces: " << m_faces.size() << std::endl;
    
    std::vector<Edges> face_edges(m_faces.size());
    for (std::size_t i = 0; i < m_faces.size(); ++i) {
      const auto& face = m_faces[i];
      auto& edges = face_edges[i];
      create_face_edges(face, edges);
    }
    update_edge_neighbors(face_edges);
    
    Edges graph_edges;
    create_graph_edges(face_edges, graph_edges);
    compute_edge_weights(graph_edges);

    create_root_level();

    if (m_faces.size() > 3)
      create_base_level_with_graphcut(
        graph_edges, face_edges);
    else create_base_level_naive(face_edges);
  }

  void compute_face_weights() {

    /* compute_naive_face_weights(); */
    compute_sum_normalized_face_weights();
    /* compute_max_normalized_face_weights(); */

    /*
    for (auto& face : m_faces)
      std::cout << "fw: " << face.weight << std::endl; */
  }

  void compute_naive_face_weights() {
    for (auto& face : m_faces)
      face.weight = 1.0;
  }

  void compute_sum_normalized_face_weights() {

    double sum = 0.0;
    for (auto& face : m_faces) {
      face.weight = get_face_weight(face);
      sum += face.weight;
    }

    if (sum != 0.0) {
      for (auto& face : m_faces)
        face.weight /= sum;
    }
  }

  void compute_max_normalized_face_weights() {

    double maxv = -1.0;
    for (auto& face : m_faces) {
      face.weight = get_face_weight(face);
      maxv = CGAL::max(maxv, face.weight);
    }

    if (maxv != 0.0) {
      for (auto& face : m_faces)
        face.weight /= maxv;
    }
  }

  double get_face_weight(const Face& face) {
    return face.get_area();
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

  void merge_edges(
    const std::vector<Edges>& face_edges,
    Edges& all_edges) {

    all_edges.clear(); std::size_t numa = 0, numb = 0;
    for (std::size_t i = 0; i < face_edges.size(); ++i) {
      for (std::size_t j = 0; j < face_edges[i].size(); ++j) {

        ++numa;
        const auto& edge = face_edges[i][j];
        const std::size_t f1 = edge.faces.first;
        const std::size_t f2 = edge.faces.second;
        if (f1 == std::size_t(-1) || f2 == std::size_t(-1)) {
          ++numb; continue;
        }

        if (f1 < f2)
          all_edges.push_back(edge);
      }
    }

    std::cout << "edges: " << 
      numa - numb << " = " << all_edges.size() << std::endl;
  }

  void create_graph_edges(
    const std::vector<Edges>& face_edges,
    Edges& graph_edges) {

    Edge edge;
    graph_edges.clear();
    for (const auto& face : m_faces) {
      edge.faces.first = face.index;
      for (const std::size_t neighbor : face.neighbors) {
        edge.faces.second = neighbor;
        if (edge.faces.first < edge.faces.second)
          graph_edges.push_back(edge);
      }
    }

    for (auto& edge : graph_edges) {
      const std::size_t f1 = edge.faces.first;
      const std::size_t f2 = edge.faces.second;

      edge.length = 
        compute_edge_length(f1, f2, face_edges);
    }

    std::sort(graph_edges.begin(), graph_edges.end(), 
    [](const Edge& a, const Edge& b) -> bool { 
      return a.get_length() > b.get_length();
    });

    std::cout << "edges: " << graph_edges.size() << std::endl;
    /*
    for (const auto& edge : gc_edges)
      std::cout << "tl: " << edge.get_length() << std::endl; */
  }

  FT compute_edge_length(
    const std::size_t f1, const std::size_t f2,
    const std::vector<Edges>& face_edges) {

    FT total_length = FT(0);
    const auto& edges = face_edges[f1];
    for (const auto& edge : edges) {
      if (edge.faces.first == f1 && edge.faces.second == f2)
        total_length += edge.compute_length();
      if (edge.faces.first == f2 && edge.faces.second == f1)
        total_length += edge.compute_length();
    }
    return total_length;
  }

  void compute_edge_weights(
    Edges& graph_edges) {

    /* compute_naive_edge_weights(graph_edges); */
    compute_sum_normalized_edge_weights(graph_edges);
    /* compute_max_normalized_edge_weights(graph_edges); */

    /*
    for (const auto& edge : graph_edges)
      std::cout << "ew: " << edge.weight << std::endl; */
  }

  void compute_naive_edge_weights(
    Edges& graph_edges) {

    for (auto& edge : graph_edges)
      edge.weight = 1.0;
  }

  void compute_sum_normalized_edge_weights(
    Edges& graph_edges) {

    double sum = 0.0;
    for (auto& edge : graph_edges) {
      edge.weight = get_edge_weight(edge); 
      if (edge.weight >= 0.0)
        sum += edge.weight;
    }

    for (auto& edge : graph_edges) {
      if (edge.weight < 0.0) {
        edge.weight = 0.0;
      } else {
        if (sum != 0.0)
          edge.weight /= sum;
      }
    }
  }

  void compute_max_normalized_edge_weights(
    Edges& graph_edges) {

    double maxv = -1.0;
    for (auto& edge : graph_edges) {
      edge.weight = get_edge_weight(edge); 
      if (edge.weight >= 0.0)
        maxv = CGAL::max(maxv, edge.weight);
    }

    for (auto& edge : graph_edges) {
      if (edge.weight < 0.0) {
        edge.weight = 0.0;
      } else {
        if (maxv != 0.0)
          edge.weight /= maxv;
      }
    }
  }

  double get_edge_weight(const Edge& edge) {
    return get_edge_weight_simple(edge);
  }

  double get_edge_weight_simple(const Edge& edge) {
    return edge.get_length();
  }

  double get_edge_weight_complex(const Edge& edge) {

    if (edge.type == Edge_type::BOUNDARY)
      return edge.get_length();

    const auto& s = edge.segment.source();
    const auto& t = edge.segment.target();

    const auto mid = internal::middle_point_2(s, t);

    const std::size_t f1 = edge.faces.first;
    const std::size_t f2 = edge.faces.second;

    const auto& plane1 = m_plane_map.at(m_faces[f1].label);
    const auto& plane2 = m_plane_map.at(m_faces[f2].label);

    const FT z1 = internal::position_on_plane_3(mid, plane1).z();
    const FT z2 = internal::position_on_plane_3(mid, plane2).z();

    const FT diff = CGAL::abs(z1 - z2);
    if (diff < m_max_height_difference / FT(2))
      return edge.get_length();
    else
      return -1.0;
  }

  void create_root_level() {
    Indices root_indices(1, 0);
    m_tree.levels.push_back(root_indices);
  }

  void create_base_level_with_graphcut(
    const Edges& graph_edges,
    const std::vector<Edges>& face_edges) {

    create_label_mappings();

    Indices gc_labels;
    set_initial_labels(gc_labels);

    std::vector<Size_pair> gc_edges;
    std::vector<double> gc_edge_weights;
    set_graphcut_edges(
      graph_edges, gc_edges, gc_edge_weights);

    std::vector< std::vector<double> > gc_cost_matrix;
    set_cost_matrix(face_edges, gc_cost_matrix);

    compute_graphcut(
      gc_edges, gc_edge_weights, gc_cost_matrix, gc_labels);
    apply_new_labels(gc_labels);
  }

  void create_label_mappings() {
    
    m_dr_mapping.clear();
    m_op_mapping.clear();

    std::size_t count = 0;
    for (const auto& face : m_faces) {
      if (m_dr_mapping.find(face.label) == m_dr_mapping.end()) {
        m_dr_mapping[face.label] = count;
        m_op_mapping[count] = face.label; 
        ++count;
      }
    }
  }

  void set_initial_labels(
    Indices& gc_labels) {

    gc_labels.clear();
    gc_labels.resize(m_faces.size());
    for (std::size_t i = 0; i < m_faces.size(); ++i) {
      gc_labels[i] = m_dr_mapping.at(m_faces[i].label);
      /* std::cout << gc_labels[i] << std::endl; */
    }
  }

  void set_graphcut_edges(
    const Edges& graph_edges,
    std::vector<Size_pair>& gc_edges,
    std::vector<double>& gc_edge_weights) {

    gc_edges.clear();
    gc_edges.reserve(graph_edges.size());

    gc_edge_weights.clear();
    gc_edge_weights.reserve(graph_edges.size());

    for (const auto& edge : graph_edges) {
      gc_edges.push_back(edge.faces);
      gc_edge_weights.push_back(get_edge_cost(edge));

      /*
      std::cout << 
      int(edge.type) << " : " <<
      gc_edges.back().first << " " << gc_edges.back().second << " " << 
      gc_edge_weights.back() << std::endl; */
    }
  }

  void set_cost_matrix(
    const std::vector<Edges>& face_edges,
    std::vector< std::vector<double> >& gc_cost_matrix) {

    const std::size_t num_labels = m_dr_mapping.size();

    gc_cost_matrix.clear();
    gc_cost_matrix.resize(num_labels);

    for (std::size_t i = 0; i < num_labels; ++i)
      gc_cost_matrix[i].resize(m_faces.size());

    std::vector<double> probabilities;
    for (std::size_t i = 0; i < m_faces.size(); ++i) {
      const auto& face = m_faces[i];
      
      probabilities.clear();
      probabilities.resize(num_labels, 0.0);

      create_probabilities(
        face, face_edges[i], probabilities);
      for (std::size_t j = 0; j < num_labels; ++j)
        gc_cost_matrix[j][i] = get_face_cost(face, probabilities[j]);
    }
  }

  double get_edge_cost(const Edge& edge) {
    return CGAL::to_double(m_beta) * edge.weight;
  }

  void create_probabilities(
    const Face& face,
    const Edges& edges,
    std::vector<double>& probabilities) {

    /*
    create_first_probabilities(0, face, edges, probabilities);
    create_naive_probabilities(face, edges, probabilities);
    create_neighbor_probabilities(face, edges, probabilities);
    create_length_probabilities(face, edges, probabilities);
    create_big_length_probabilities(face, edges, probabilities); */

    create_stable_probabilities(face, edges, probabilities);
    normalize_with_sum(probabilities);

    /*
    for (const auto& value : probabilities)
      std::cout << value << " ";
    std::cout << std::endl; */
  }

  void normalize_with_sum(
    std::vector<double>& probabilities) {
    
    double sum = 0.0;
    for (const auto& value : probabilities)
      sum += value;

    if (sum != 0.0) {
      for (auto& value : probabilities)
        value /= sum;
    }
  }

  void normalize_with_max(
    std::vector<double>& probabilities) {
    
    double maxv = -1.0;
    for (const auto& value : probabilities)
      maxv = CGAL::max(maxv, value);

    if (maxv != 0.0) {
      for (auto& value : probabilities)
        value /= maxv;
    }
  }

  void create_first_probabilities(
    const std::size_t fi,
    const Face& face,
    const Edges& edges,
    std::vector<double>& probabilities) {

    const std::size_t idx = m_dr_mapping.at(m_faces[fi].label);
    probabilities[idx] = 1.0;
  }

  void create_naive_probabilities(
    const Face& face,
    const Edges& edges,
    std::vector<double>& probabilities) {

    const std::size_t idx = m_dr_mapping.at(face.label);
    probabilities[idx] = 1.0;
  }

  void create_neighbor_probabilities(
    const Face& face,
    const Edges& edges,
    std::vector<double>& probabilities) {

    for (const auto& edge : edges) {
      if (edge.type == Edge_type::BOUNDARY)
        continue;

      const std::size_t nf = edge.faces.second;
      const auto& nface = m_faces[nf];
      const std::size_t idx = m_dr_mapping.at(nface.label);
      probabilities[idx] += 1.0;
    } 
  }

  void create_length_probabilities(
    const Face& face,
    const Edges& edges,
    std::vector<double>& probabilities) {

    for (const auto& edge : edges) {
      if (edge.type == Edge_type::BOUNDARY)
        continue;

      const std::size_t nf = edge.faces.second;
      const auto& nface = m_faces[nf];
      const std::size_t idx = m_dr_mapping.at(nface.label);
      const double length = static_cast<double>(edge.compute_length());
      probabilities[idx] += 1.0 * length;
    } 
  }

  void create_big_length_probabilities(
    const Face& face,
    const Edges& edges,
    std::vector<double>& probabilities) {

    for (const auto& edge : edges) {
      if (edge.type == Edge_type::BOUNDARY)
        continue;

      const std::size_t f1 = edge.faces.first;
      const std::size_t f2 = edge.faces.second;

      const auto& face1 = m_faces[f1];
      const auto& face2 = m_faces[f2];

      const std::size_t idx1 = m_dr_mapping.at(face1.label);
      const std::size_t idx2 = m_dr_mapping.at(face2.label);

      const double length = static_cast<double>(edge.compute_length());

      if (face1.get_area() > face2.get_area())
        probabilities[idx1] += 1.0 * length;
      else
        probabilities[idx2] += 1.0 * length;
    } 
  }

  void create_stable_probabilities(
    const Face& face,
    const Edges& edges,
    std::vector<double>& probabilities) {

    double length = 0.0;
    for (const auto& edge : edges)
      if (comply_by_angle(edge))
        length += edge.compute_length();
    if (length >= CGAL::to_double(m_min_length_2)) {
      probabilities[m_dr_mapping.at(face.label)] = 1.0; return;
    }

    for (const auto& edge : edges) {
      if (edge.type != Edge_type::BOUNDARY)
        add_internal_edge_weight(edge, probabilities);
    } 
  }

  void add_boundary_edge_weight(
    const Edge& edge,
    std::vector<double>& probabilities) {

    const double length = 
      static_cast<double>(edge.compute_length());
    const std::size_t i1 = get_first_index(edge);

    probabilities[i1] += length * get_basic_rate(edge);
    probabilities[i1] += length * get_corner_rate(edge);
    probabilities[i1] += length * get_parallelism_rate(edge);
  }

  void add_internal_edge_weight(
    const Edge& edge,
    std::vector<double>& probabilities) {

    const double length = 
      static_cast<double>(edge.compute_length());

    const std::size_t i1 = get_first_index(edge);
    const std::size_t i2 = get_second_index(edge);
    const std::size_t ib = get_biggest_face_index(edge);

    probabilities[ib] += length * get_basic_rate(edge);

    /*
    probabilities[ib] += length * get_adjacency_rate(edge);
    probabilities[ib] += length * get_coplanarity_rate(edge);
    probabilities[i1] += length * get_parallelism_rate(edge); */
  }

  std::size_t get_first_index(const Edge& edge) {

    const std::size_t f1 = edge.faces.first;
    const auto& face1 = m_faces[f1];
    return m_dr_mapping.at(face1.label);
  }

  std::size_t get_second_index(const Edge& edge) {

    const std::size_t f2 = edge.faces.second;
    const auto& face2 = m_faces[f2];
    return m_dr_mapping.at(face2.label);
  }

  std::size_t get_biggest_face_index(const Edge& edge) {

    const std::size_t f1 = edge.faces.first;
    const std::size_t f2 = edge.faces.second;

    const auto& face1 = m_faces[f1];
    const auto& face2 = m_faces[f2];

    const std::size_t idx1 = m_dr_mapping.at(face1.label);
    const std::size_t idx2 = m_dr_mapping.at(face2.label);

    if (face1.get_area() > face2.get_area()) return idx1;
    else return idx2;
  }

  double get_basic_rate(const Edge& edge) {
    return 1.0;
  }

  double get_corner_rate(const Edge& edge) {

    const auto& v1 = m_vertices[edge.from_vertex];
    const auto& v2 = m_vertices[edge.to_vertex];

    if (
      v1.type == Point_type::OUTER_CORNER && 
      v2.type == Point_type::OUTER_CORNER) return 1.0;
    return 0.0;
  }

  double get_adjacency_rate(const Edge& edge) {
    
    const std::size_t f1 = edge.faces.first;
    const std::size_t f2 = edge.faces.second;

    const auto& face1 = m_faces[f1];
    const auto& face2 = m_faces[f2];

    const auto& plane1 = m_plane_map.at(face1.label);
    const auto& plane2 = m_plane_map.at(face2.label);

    const auto& s = edge.segment.source();
    const auto& t = edge.segment.target();

    const auto mid = internal::middle_point_2(s, t);
    
    const FT z1 = internal::position_on_plane_3(mid, plane1).z();
    const FT z2 = internal::position_on_plane_3(mid, plane2).z();

    const FT diff = CGAL::abs(z1 - z2);
    if (diff < m_max_height_difference / FT(2)) return 1.0;
    return 0.0;
  }

  double get_coplanarity_rate(const Edge& edge) {

    const std::size_t f1 = edge.faces.first;
    const std::size_t f2 = edge.faces.second;

    const auto& face1 = m_faces[f1];
    const auto& face2 = m_faces[f2];

    const auto& plane1 = m_plane_map.at(face1.label);
    const auto& plane2 = m_plane_map.at(face2.label);

    auto normal1 = plane1.orthogonal_vector();
    internal::normalize(normal1);
    auto normal2 = plane2.orthogonal_vector();
    internal::normalize(normal2);

    const FT angle_deg = 
      CGAL::abs(internal::angle_3d(normal1, normal2));
    if (CGAL::abs(angle_deg) < m_bound_min) return 1.0;
    return 0.0;
  }

  double get_parallelism_rate(const Edge& edge) {

    if (comply_by_angle_and_distance(edge)) 
      return 1.0;
    return 0.0;
  }

  bool comply_by_angle(
    const Edge& edge) {

    for (const auto& direction : m_directions) {
      const FT angle = angle_degree_2(
        direction, edge.segment);
      const FT angle_2 = get_angle_2(angle);

      if ( 
        (CGAL::abs(angle_2) <= m_bound_min) ||
        (CGAL::abs(angle_2) >= m_bound_max) ) 
        return true;
    }
    return false;
  }

  bool comply_by_angle_and_distance(
    const Edge& edge) {

    const FT elength = edge.compute_length();
    for (const auto& direction : m_directions) {   
      const FT dlength = 
        internal::distance(direction.source(), direction.target());
      if (elength < dlength / FT(2))
        continue;

      const FT angle = angle_degree_2(
        direction, edge.segment);
      const FT angle_2 = get_angle_2(angle);

      if ( 
        (CGAL::abs(angle_2) <= m_bound_min) ||
        (CGAL::abs(angle_2) >= m_bound_max) ) 
        return true;
    }
    return false;
  }

  double get_face_cost(
    const Face& face,
    const double probability) {
    
    return (1.0 - probability) * face.weight;
  }

  void compute_graphcut(
    const std::vector<Size_pair>& gc_edges,
    const std::vector<double>& gc_edge_weights,
    const std::vector< std::vector<double> >& gc_cost_matrix,
    std::vector<std::size_t>& gc_labels) {

    Alpha_expansion graphcut;
    graphcut(
      gc_edges, gc_edge_weights, gc_cost_matrix, gc_labels);
    std::cout << "gc finished" << std::endl;
  }

  void apply_new_labels(
    const Indices& gc_labels) {

    for (std::size_t i = 0; i < gc_labels.size(); ++i) {
      auto& face = m_faces[i];
      face.label = m_op_mapping.at(gc_labels[i]);
      for (auto fh = face.tri.delaunay.finite_faces_begin();
      fh != face.tri.delaunay.finite_faces_end(); ++fh)
        fh->info().label = face.label;
    }
  }

  void create_base_level_naive(
    const std::vector<Edges>& face_edges) {

    auto& nodes  = m_tree.nodes;
    auto& levels = m_tree.levels;
    
    const FT avg_area = compute_average_face_area();
    const FT eps = avg_area / FT(4);

    Indices level;
    for (std::size_t i = 0; i < face_edges.size(); ++i) {
      const auto& edges = face_edges[i];
      auto& face = m_faces[i];
      if (face.area < eps) continue;

      for (const auto& edge : edges) {
        if (comply_by_angle_and_distance(edge)) {
          level.push_back(i + 1); 
          face.level = 1; break;
        }
      }
    }

    levels.push_back(level);
    auto& root = nodes[0];
    root.children = level;
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

  void create_sublevels_naive(
    const std::vector<Edges>& face_edges) {
    
    Indices level;
    auto& levels = m_tree.levels;

    bool completed = false;
    std::size_t deep = 1;

    do {
      deep += 1;
      completed = create_sublevel_naive(
        deep, face_edges, level);
      levels.push_back(level);
    } while (!completed && deep <= 100);

    if (deep > 100) {
      std::cerr << 
      "ERROR: too many iterations, create_sublevels()!" << std::endl;
      exit(1);
    }
  }

  bool create_sublevel_naive(
    const std::size_t deep,
    const std::vector<Edges>& face_edges,
    Indices& level) {

    level.clear();
    bool completed = true;
    auto& nodes = m_tree.nodes;

    for (std::size_t i = 0; i < face_edges.size(); ++i) {
      const auto& edges = face_edges[i];
      auto& face = m_faces[i];
      if (face.level != std::size_t(-1)) continue;
      
      completed = false;
      const auto& neighbors = face.neighbors;
      const std::size_t best_face = get_best_face(
        deep, edges, face.label, neighbors);

      if (best_face == std::size_t(-1)) {
        /* auto& base = m_tree.levels[1];
        base.push_back(i + 1);
        auto& root = nodes[0];
        root.children.push_back(i + 1);
        face.level = 1; */
        continue;
      }

      const std::size_t node_idx = best_face + 1;
      level.push_back(i + 1);
      nodes[node_idx].children.push_back(i + 1);
      face.level = deep;
    }
    return completed;
  }

  std::size_t get_best_face(
    const std::size_t deep,
    const Edges& edges,
    const std::size_t query_label,
    const Indices& neighbors) {

    if (query_label == std::size_t(-1)) 
      return std::size_t(-1);
    const auto& nodes = m_tree.nodes;

    std::map<std::size_t, FT> length;
    for (const std::size_t idx : neighbors) {
      if (m_faces[idx].level != deep - 1) continue;
      length[idx] = FT(0);
    }

    for (const auto& edge : edges) {
      const std::size_t f2 = edge.faces.second;
      const auto& seg = edge.segment;

      if (length.find(f2) != length.end())
        length[f2] += internal::distance(seg.source(), seg.target());
    }

    const auto& query_plane = m_plane_map.at(query_label);
    auto query_normal = query_plane.orthogonal_vector();
    internal::normalize(query_normal);

    std::size_t best_idx = std::size_t(-1);
    FT min_angle = internal::max_value<FT>();
    for (const auto& pair : length) {
      if (pair.second > FT(0)) {

        const std::size_t idx = pair.first;
        const auto& node = nodes[idx + 1];
        const auto& plane = m_plane_map.at(node.label);
        auto normal = plane.orthogonal_vector();
        internal::normalize(normal);

        const FT angle_deg = 
          CGAL::abs(internal::angle_3d(query_normal, normal));
        if (angle_deg < min_angle) {
          min_angle = angle_deg;
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

    for (auto& face : m_faces) {
      for (auto fh = face.tri.delaunay.finite_faces_begin(); 
      fh != face.tri.delaunay.finite_faces_end(); ++fh)
        fh->info().label = face.label;
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
