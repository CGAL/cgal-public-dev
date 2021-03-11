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
#include <CGAL/Levels_of_detail/internal/Shape_detection/Region_growing.h>

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
typename Edge_type,
typename Face_type>
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
  using Triangle_2 = typename Traits::Triangle_2;

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
  using Triangulation = internal::Triangulation<Traits>;
  using LF_circulator = typename Triangulation::Delaunay::Line_face_circulator;

  class Face_neighbor_query {

  public:
    Face_neighbor_query(
      const std::vector<Face>& faces,
      const Indices& indices) :
    m_faces(faces),
    m_indices(indices)
    { }

    void operator()(
      const std::size_t query_index,
      std::vector<std::size_t>& neighbors) const {

      CGAL_assertion(query_index >= 0);
      CGAL_assertion(query_index < m_indices.size());
      neighbors.clear();
      const auto& ns = m_faces[m_indices[query_index]].neighbors;
      for (const std::size_t idx : ns) {
        for (std::size_t i = 0; i < m_indices.size(); ++i) {
          if (m_indices[i] == idx) {
            neighbors.push_back(i); break;
          }
        }
      }
    }

  private:
    const std::vector<Face>& m_faces;
    const Indices& m_indices;
  };

  class Face_region {
  public:
    Face_region(
      const std::vector<Face>& faces,
      const Indices& indices,
      const std::size_t min_region_size = 1) :
    m_faces(faces),
    m_indices(indices),
    m_min_region_size(min_region_size)
    { }

    bool is_already_visited(
      const std::size_t,
      const std::size_t query_index,
      const bool is_visited) const { return false; }

    bool is_part_of_region(
      const std::size_t, const std::size_t,
      const std::vector<std::size_t>&) const {
      return true;
    }

    bool is_valid_region(const std::vector<std::size_t>& region) const {
      return ( region.size() >= m_min_region_size );
    }

    void update(const std::vector<std::size_t>&) { }

  private:
    const std::vector<Face>& m_faces;
    const Indices& m_indices;
    const std::size_t m_min_region_size;
  };

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
        if (child == std::size_t(-1)) continue;
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
      if (node_idx == std::size_t(-1)) return;
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
    const FT beta,
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
  m_beta(beta) {

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
    if (m_faces.size() <= 3)
      return;
    set_original_labels();
    build_tree_with_graphcut();
    set_tree_leaves();
  }

  void build_v1() {
    set_original_labels();
    build_tree_naive();
    set_tree_leaves();
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
    if (m_faces.size() <= 3)
      return;
    if (level < 0) level = 0;
    if (level >= m_tree.levels.size())
      level = m_tree.levels.size() - 1;
    cut_along_tree(level);
  }

  void merge_faces() {

    update_edge_neighbors();
    std::set<std::size_t> labels;
    for (const auto& face : m_faces)
      if (face.label != std::size_t(-1))
        labels.insert(face.label);

    Face face;
    std::vector<Face> faces;
    std::size_t count = 0;

    if (labels.size() == 1) {

      const std::size_t ref_label = *(labels.begin());
      face.label = ref_label;
      traverse(ref_label, false, face.hedges);
      const bool success = initialize_face(face);
      if (success) {
        face.index = count; ++count;
        faces.push_back(face);
      }

    } else {

      for (const std::size_t ref_label : labels) {
        bool finished = false;
        set_default_edges();

        do {
          face.label = ref_label;
          finished = traverse(ref_label, true, face.hedges);
          if (!finished) {
            const bool success = initialize_face(face);
            if (success) {
              face.index = count; ++count;
              faces.push_back(face);
            }
          }
        } while (!finished);
        set_default_edges();
      }
    }

    m_faces = faces;
    sort_faces();
    create_face_neighbors();
    set_face_types();
    update_edge_neighbors();
    update_corners();

    std::cout << "num merged faces: " << m_faces.size() << std::endl;
  }

  void remove_one_neighbor_faces() {

    if (m_faces.size() <= 3)
      return;

    const FT avg_area = compute_average_face_area();
    const FT eps = avg_area / FT(4);

    for (auto& face : m_faces) {
      if (face.neighbors.size() == 1 && face.area < eps) {
        face.label = m_faces[face.neighbors[0]].label;
        for (auto fh = face.tri.delaunay.finite_faces_begin();
        fh != face.tri.delaunay.finite_faces_end(); ++fh)
          fh->info().label = face.label;
      }
    }
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

  void update_corners() {

    std::set<std::size_t> unique;
    for (auto& vertex : m_vertices) {
      if (vertex.type == Point_type::CORNER) {

        unique.clear();
        for (const std::size_t he_idx : vertex.hedges) {
          const auto& edge = m_edges[m_halfedges[he_idx].edg_idx];

          const auto& faces = edge.faces;
          const std::size_t f1 = faces.first;
          const std::size_t f2 = faces.second;
          if (f1 == std::size_t(-1) || f2 == std::size_t(-1))
            continue;

          const std::size_t l1 = m_faces[f1].label;
          const std::size_t l2 = m_faces[f2].label;
          if (l1 == std::size_t(-1) || l2 == std::size_t(-1))
            continue;

          if (l1 != l2) {
            unique.insert(l1);
            unique.insert(l2);
          }
        }

        if (unique.size() < 3)
          vertex.type = Point_type::LINEAR;
      }
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

  bool traverse(
    const std::size_t ref_label,
    const bool skip_bounds,
    Indices& hedges) {

    hedges.clear();
    for (const auto& he : m_halfedges) {
      const auto& edge = m_edges[he.edg_idx];
      if (edge.used) continue;

      const std::size_t f1 = edge.faces.first;
      const std::size_t f2 = edge.faces.second;

      if (skip_bounds) {
        if (f1 == std::size_t(-1) || f2 == std::size_t(-1))
          continue;
      }

      if (f1 == std::size_t(-1) && f2 != std::size_t(-1)) {
        traverse(he.index, ref_label, hedges);
        return false;
      }

      if (f2 == std::size_t(-1) && f1 != std::size_t(-1)) {
        traverse(he.index, ref_label, hedges);
        return false;
      }

      const auto& face1 = m_faces[f1];
      const auto& face2 = m_faces[f2];

      if (face1.label != ref_label && face2.label != ref_label)
        continue;

      if (face1.label == ref_label && face2.label != ref_label) {
        traverse(he.index, ref_label, hedges);
        return false;
      }
      if (face2.label == ref_label && face1.label != ref_label) {
        traverse(he.index, ref_label, hedges);
        return false;
      }
    }
    return true;
  }

  void traverse(
    const std::size_t idx,
    const std::size_t ref_label,
    Indices& hedges) {

    set_default_next_halfedges();
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
      m_edges[m_halfedges[curr].edg_idx].used = true;

      auto& other = m_halfedges[curr];
      const std::size_t to_idx = other.to_vertex;
      const auto& to = m_vertices[to_idx];
      find_next(start, ref_label, to, other);

      const auto& s = m_vertices[other.from_vertex].point;
      const auto& t = m_vertices[other.to_vertex].point;

      segments.push_back(Segment_2(s, t));
      curr = other.next;

      if (count >= 10000) {

        std::cout.precision(30);
        std::cout << "Error: traverse() max count reached!" << std::endl;
        std::cout << "Ref label: " << ref_label << std::endl;
        for (const auto& face : m_faces)
          std::cout << face.label << " " << face.skip << " " << face.area << std::endl;

        Saver saver;
        saver.save_polylines(
        segments, "/Users/monet/Documents/gf/lod/logs/buildings/tmp/debug-count");
        exit(EXIT_FAILURE);
      }

      if (curr == std::size_t(-1)) {

        std::cout.precision(30);
        std::cout << "Error: traverse() failed!" << std::endl;
        std::cout << "Ref label: " << ref_label << std::endl;
        for (const auto& face : m_faces)
          std::cout << face.label << " " << face.skip << " " << face.area << std::endl;

        Saver saver;
        saver.save_polylines(
        segments, "/Users/monet/Documents/gf/lod/logs/buildings/tmp/debug-fail");
        exit(EXIT_FAILURE);
      }
      ++count;

    } while (curr != start);
    set_default_next_halfedges();
  }

  void set_default_next_halfedges() {
    for (auto& he : m_halfedges)
      he.next = std::size_t(-1);
  }

  void set_default_edges() {
    for (auto& edge : m_edges)
      edge.used = false;
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
      const auto& faces = edge.faces;
      const std::size_t f1 = faces.first;
      const std::size_t f2 = faces.second;

      std::size_t l1 = std::size_t(-1);
      if (f1 != std::size_t(-1))
        l1 = m_faces[f1].label;

      std::size_t l2 = std::size_t(-1);
      if (f2 != std::size_t(-1))
        l2 = m_faces[f2].label;

      if (l1 == ref_label && l2 != ref_label) {
        he.next = other.index;
        return;
      }
      if (l2 == ref_label && l1 != ref_label) {
        he.next = other.index;
        return;
      }
    }
  }

  bool initialize_face(
    Face& face) {

    remove_collinear_edges(face);
    const bool success = create_face_triangulation(face);
    if (!success) return false;
    create_face_visibility(face);
    update_face_label(face);
    compute_face_area(face);
    return true;
  }

  void remove_collinear_edges(Face& face) {

    const std::size_t m = face.hedges.size();
    for (std::size_t i = 0; i < m; ++i) {

      const std::size_t hi  = face.hedges[i];
      auto& curr = m_vertices[m_halfedges[hi].from_vertex];
      const auto& edge = m_edges[m_halfedges[hi].edg_idx];

      const std::size_t f1 = edge.faces.first;
      const std::size_t f2 = edge.faces.second;

      if (f1 != std::size_t(-1) && f2 != std::size_t(-1))
        continue;
      if (curr.type == Point_type::OUTER_CORNER)
        continue;

      bool skip = true;
      for (const std::size_t he_idx : curr.hedges) {
        const auto& tmp = m_edges[m_halfedges[he_idx].edg_idx];

        const std::size_t ff1 = tmp.faces.first;
        const std::size_t ff2 = tmp.faces.second;

        if (ff1 != std::size_t(-1)) {
          const auto& fface1 = m_faces[ff1];
          if (fface1.label != face.label) {
            skip = false; break;
          }
        }

        if (ff2 != std::size_t(-1)) {
          const auto& fface2 = m_faces[ff2];
          if (fface2.label != face.label) {
            skip = false; break;
          }
        }
      }

      curr.skip = skip;
    }
  }

  bool create_face_triangulation(Face& face) {

    auto& tri = face.tri.delaunay;
    tri.clear();

    std::vector<Edge> edges;
    create_face_edges(face, edges);
    if (edges.size() == 0) {
      std::cout <<
      "Warning: empty face, create_face_triangulation(), image tree!"
      << std::endl;
      return false;
    }

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

  void set_original_labels() {
    for (auto& face : m_faces)
      face.original = face.label;
  }

  void set_tree_leaves() {
    for (auto& node : m_tree.nodes)
      if (node.children.size() == 0)
        node.type = Node_type::LEAF;
  }

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
    segments, "/Users/monet/Documents/gf/lod/logs/buildings/tmp/skeleton");
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
    create_base_level_naive(
      face_edges);
    create_sublevels_naive(
      face_edges);
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
    create_base_level_with_graphcut(
      graph_edges, face_edges);
    create_sublevels_with_graphcut(
      face_edges);
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

  void update_edge_neighbors() {

    for (auto& edge : m_edges) {
      edge.faces.first  = std::size_t(-1);
      edge.faces.second = std::size_t(-1);
    }

    for (std::size_t i = 0; i < m_faces.size(); ++i) {
      const auto& face = m_faces[i];
      const auto& hedges = face.hedges;

      for (const std::size_t he_idx : hedges) {
        const auto& he = m_halfedges[he_idx];
        auto& edge = m_edges[he.edg_idx];

        const std::size_t v1 = edge.from_vertex;
        const std::size_t v2 = edge.to_vertex;

        edge.faces.first  = i;
        edge.faces.second = find_neighbor_face(i, v1, v2);
        if (edge.faces.second == std::size_t(-1))
          edge.type = Edge_type::BOUNDARY;
        else
          edge.type = Edge_type::INTERNAL;
      }
    }
  }

  std::size_t find_neighbor_face(
    const std::size_t skip,
    const std::size_t v1, const std::size_t v2) {

    for (std::size_t i = 0; i < m_faces.size(); ++i) {
      if (i == skip) continue;

      const auto& face = m_faces[i];
      const auto& hedges = face.hedges;

      for (const std::size_t he_idx : hedges) {
        const auto& he = m_halfedges[he_idx];
        const auto& edge = m_edges[he.edg_idx];

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
    add_base_level_from_graphcut();
  }

  void add_base_level_from_graphcut() {

    auto& nodes  = m_tree.nodes;
    auto& levels = m_tree.levels;

    std::set<std::size_t> labels;
    for (const auto& face : m_faces)
      labels.insert(face.label);

    Indices indices, level;
    for (const std::size_t label : labels) {
      indices.clear();
      for (std::size_t i = 0; i < m_faces.size(); ++i)
        if (m_faces[i].label == label)
          indices.push_back(i);

      std::sort(indices.begin(), indices.end(),
      [&](const std::size_t id1, const std::size_t id2) -> bool {
        return m_faces[id1].area > m_faces[id2].area;
      });

      for (const std::size_t idx : indices) {
        if (m_faces[idx].label == m_faces[idx].original) {
          level.push_back(idx + 1); m_faces[idx].level = 1;
          break;
        }
      }
    }

    levels.push_back(level);
    auto& root = nodes[0];
    root.children = level;
  }

  void create_sublevels_with_graphcut(
    const std::vector<Edges>& face_edges) {

    add_mansards(face_edges);
    add_mansard_sublevels(face_edges);
  }

  void add_mansards(
    const std::vector<Edges>& face_edges) {

    using Region_growing = internal::Region_growing<
    Indices, Face_neighbor_query, Face_region>;

    Indices indices;
    std::vector<Indices> regions;

    auto& nodes  = m_tree.nodes;
    auto& levels = m_tree.levels;

    Indices level;
    std::size_t count = nodes.size() - 1;

    const auto& base = levels[1];
    for (const std::size_t idx : base) {
      const std::size_t ref_idx = idx - 1;

      indices.clear();
      for (const auto& face : m_faces) {
        if (face.index == ref_idx) continue;
        if (face.label == m_faces[ref_idx].label)
          indices.push_back(face.index);
      }

      if (indices.size() == 0) continue;
      Face_neighbor_query neighbor_query(
        m_faces, indices);
      Face_region region_type(
        m_faces, indices);
      Region_growing region_growing(
        indices, neighbor_query, region_type);

      regions.clear();
      region_growing.detect(std::back_inserter(regions));

      Node node;
      for (const auto& region : regions) {
        node.children.clear();
        for (const std::size_t ri : region)
          node.children.push_back(indices[ri] + 1);

        node.index = count + 1; ++count;
        node.face_index = std::size_t(-1);
        node.label = std::size_t(-1);
        node.type = Node_type::CHILD;
        nodes.push_back(node);

        level.push_back(node.index);
        nodes[idx].children.push_back(node.index);
      }
    }
    levels.push_back(level);
  }

  void add_mansard_sublevels(
    const std::vector<Edges>& face_edges) {

    auto& nodes  = m_tree.nodes;
    auto& levels = m_tree.levels;
    const auto& base = levels[2];

    Indices level;
    for (const std::size_t idx : base) {
      for (const std::size_t child : nodes[idx].children)
        level.push_back(child);
    }
    levels.push_back(level);
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

    if (m_faces.size() == 1)
      return;

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

    if (m_faces.size() == 1)
      return;

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

    if (lidx == std::size_t(-1)) return;
    const auto& nodes = m_tree.nodes;
    const auto& level = m_tree.levels[lidx];

    if (nodes.size() == 1) {
      CGAL_assertion(m_faces.size() == nodes.size());
      m_faces[0].label = nodes[0].label; return;
    }

    for (auto& face : m_faces)
      face.label = face.original;

    Indices findices;
    for (const std::size_t nidx : level) {
      if (nidx == std::size_t(-1)) continue;

      m_tree.traverse_children(nidx, findices);
      for (const std::size_t fidx : findices) {

        if (fidx != std::size_t(-1)) {
          auto& face = m_faces[fidx];
          face.label = nodes[nidx].label;
        }
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
};

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_INTERNAL_IMAGE_TREE_H
