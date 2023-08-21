// Copyright (c) 2023 GeometryFactory.
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Ken Arroyo Ohori
//                 Guillaume Damiand

#ifndef CGAL_POLYGON_REPAIR_2_H
#define CGAL_POLYGON_REPAIR_2_H

#include <CGAL/license/Polygon_repair_2.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>

#include <CGAL/Polygon_repair_2/Multipolygon_with_holes_2.h>
#include <CGAL/Polygon_repair_2/Triangulation_face_base_with_repair_info_2.h>
#include <CGAL/Polygon_repair_2/Triangulation_with_odd_even_constraints_2.h>

namespace CGAL {

namespace Polygon_repair_2 {

template <class Kernel, class PolygonContainer>
class Polygon_repair_2;

/// \ingroup PkgPolygonRepair2Functions
/// Repair a polygon without holes
template <class Kernel, class PolygonContainer>
Multipolygon_with_holes_2<Kernel, PolygonContainer> repair_odd_even(const Polygon_2<Kernel, PolygonContainer>& p) {
  CGAL::Polygon_repair_2::Polygon_repair_2<Kernel, PolygonContainer> pr;
  pr.add_to_triangulation_odd_even(p);
  if (pr.triangulation().number_of_faces() > 0) {
    pr.label_triangulation_odd_even();
    pr.reconstruct_multipolygon();
  } return pr.multipolygon();
}

/// \ingroup PkgPolygonRepair2Functions
/// Repair a polygon with holes
template <class Kernel, class PolygonContainer>
Multipolygon_with_holes_2<Kernel, PolygonContainer> repair_odd_even(const Polygon_with_holes_2<Kernel, PolygonContainer>& p) {
  CGAL::Polygon_repair_2::Polygon_repair_2<Kernel, PolygonContainer> pr;
  pr.add_to_triangulation_odd_even(p);
  if (pr.triangulation().number_of_faces() > 0) {
    pr.label_triangulation_odd_even();
    pr.reconstruct_multipolygon();
  } return pr.multipolygon();
}

/// \ingroup PkgPolygonRepair2Functions
/// Repair a multipolygon with holes
template <class Kernel, class PolygonContainer>
Multipolygon_with_holes_2<Kernel, PolygonContainer> repair_odd_even(const Multipolygon_with_holes_2<Kernel, PolygonContainer>& mp) {
  CGAL::Polygon_repair_2::Polygon_repair_2<Kernel, PolygonContainer> pr;
  pr.add_to_triangulation_odd_even(mp);
  if (pr.triangulation().number_of_faces() > 0) {
    pr.label_triangulation_odd_even();
    pr.reconstruct_multipolygon();
  } return pr.multipolygon();
}

/*! \ingroup PkgPolygonRepair2Ref
 *
 * The class `Polygon_repair_2` builds on a constrained
 * triangulation to remove the parts of constraints that overlap an even number of times
 */
template <class Kernel = CGAL::Exact_predicates_inexact_constructions_kernel,
          class PolygonContainer = std::vector<typename Kernel::Point_2>>
class Polygon_repair_2 {
public:
  using Vertex_base = CGAL::Triangulation_vertex_base_2<Kernel>;
  using Face_base = CGAL::Constrained_triangulation_face_base_2<Kernel>;
  using Face_base_with_repair_info = CGAL::Triangulation_face_base_with_repair_info_2<Kernel, Face_base>;
  using Triangulation_data_structure = CGAL::Triangulation_data_structure_2<Vertex_base, Face_base_with_repair_info>;
  using Tag = CGAL::Exact_predicates_tag; // assumed for now
  using Constrained_Delaunay_triangulation = CGAL::Constrained_Delaunay_triangulation_2<Kernel, Triangulation_data_structure, Tag>;
  using Triangulation = Triangulation_with_odd_even_constraints_2<Constrained_Delaunay_triangulation>;

  struct Polygon_less {
    using Polygon_2 = CGAL::Polygon_2<Kernel, PolygonContainer>;
    bool operator()(const Polygon_2& pa, const Polygon_2& pb) const {
      typename Polygon_2::Vertex_iterator va = pa.vertices_begin();
      typename Polygon_2::Vertex_iterator vb = pb.vertices_begin();
      while (va != pa.vertices_end() && vb != pb.vertices_end()) {
        if (*va != *vb) return *va < *vb;
        ++va;
        ++vb;
      }
      if (vb == pb.vertices_end()) return false;
      return true;
    }
  };

  struct Polygon_with_holes_less {
    using Polygon_with_holes_2 = CGAL::Polygon_with_holes_2<Kernel, PolygonContainer>;
    Polygon_less pl;
    bool operator()(const Polygon_with_holes_2& pa, const Polygon_with_holes_2& pb) const {
      if (pl(pa.outer_boundary(), pb.outer_boundary())) return true;
      if (pl(pb.outer_boundary(), pa.outer_boundary())) return false;
      typename Polygon_with_holes_2::Hole_const_iterator ha = pa.holes_begin();
      typename Polygon_with_holes_2::Hole_const_iterator hb = pb.holes_begin();
      while (ha != pa.holes_end() && hb != pb.holes_end()) {
        if (pl(*ha, *hb)) return true;
        if (pl(*hb, *ha)) return false;
      }
      if (hb == pb.holes_end()) return false;
      return true;
    }
  };

  /// \name Creation
  Polygon_repair_2() : number_of_polygons(0), number_of_holes(0) {}

  /// \name Modifiers
  /// @{

  // Add edges of the polygon to the triangulation
  void add_to_triangulation_odd_even(const Polygon_2<Kernel, PolygonContainer>& polygon) {
    std::unordered_set<std::pair<typename Kernel::Point_2, typename Kernel::Point_2>,
                       boost::hash<std::pair<typename Kernel::Point_2, typename Kernel::Point_2>>> unique_edges;

    // Get unique edges
    for (auto const& edge: polygon.edges()) {
      if (edge.source() == edge.target()) continue;
      std::pair<typename Kernel::Point_2, typename Kernel::Point_2> pair = (edge.source() < edge.target())?
      std::make_pair(edge.source(), edge.target()) : std::make_pair(edge.target(), edge.source());
      auto inserted = unique_edges.insert(pair);
      if (!inserted.second) unique_edges.erase(inserted.first);
    }

    // Insert vertices
    std::unordered_map<typename Kernel::Point_2, typename Triangulation::Vertex_handle> vertices;
    std::vector<std::pair<typename Triangulation::Vertex_handle, typename Triangulation::Vertex_handle>> edges_to_insert;
    edges_to_insert.reserve(unique_edges.size());
    for (auto const& edge: unique_edges) {
      typename Triangulation::Vertex_handle first_vertex, second_vertex;
      typename std::unordered_map<typename Kernel::Point_2, typename Triangulation::Vertex_handle>::const_iterator found = vertices.find(edge.first);
      if (found == vertices.end()) {
        first_vertex = t.insert(edge.first, search_start);
        vertices[edge.first] = first_vertex;
      } else {
        first_vertex = found->second;
      } search_start = first_vertex->face();
      found = vertices.find(edge.second);
      if (found == vertices.end()) {
        second_vertex = t.insert(edge.second, search_start);
        vertices[edge.second] = second_vertex;
      } else {
        second_vertex = found->second;
      } search_start = second_vertex->face();
      edges_to_insert.emplace_back(first_vertex, second_vertex);
    }

    // Insert edges
    for (auto const& edge: edges_to_insert) {
      t.odd_even_insert_constraint(edge.first, edge.second);
    }
  }

  // Add edges of the polygon to the triangulation
  void add_to_triangulation_odd_even(const Polygon_with_holes_2<Kernel, PolygonContainer>& polygon) {
    std::unordered_set<std::pair<typename Kernel::Point_2, typename Kernel::Point_2>,
                       boost::hash<std::pair<typename Kernel::Point_2, typename Kernel::Point_2>>> unique_edges;

    // Get unique edges
    for (auto const& edge: polygon.outer_boundary().edges()) {
      if (edge.source() == edge.target()) continue;
      std::pair<typename Kernel::Point_2, typename Kernel::Point_2> pair = (edge.source() < edge.target())?
      std::make_pair(edge.source(), edge.target()) : std::make_pair(edge.target(), edge.source());
      auto inserted = unique_edges.insert(pair);
      if (!inserted.second) unique_edges.erase(inserted.first);
    }
    for (auto const& hole: polygon.holes()) {
      for (auto const& edge: hole.edges()) {
        if (edge.source() == edge.target()) continue;
        std::pair<typename Kernel::Point_2, typename Kernel::Point_2> pair = (edge.source() < edge.target())?
        std::make_pair(edge.source(), edge.target()) : std::make_pair(edge.target(), edge.source());
        auto inserted = unique_edges.insert(pair);
        if (!inserted.second) unique_edges.erase(inserted.first);
      }
    }

    // Insert vertices
    std::unordered_map<typename Kernel::Point_2, typename Triangulation::Vertex_handle> vertices;
    std::vector<std::pair<typename Triangulation::Vertex_handle, typename Triangulation::Vertex_handle>> edges_to_insert;
    edges_to_insert.reserve(unique_edges.size());
    for (auto const& edge: unique_edges) {
      typename Triangulation::Vertex_handle first_vertex, second_vertex;
      typename std::unordered_map<typename Kernel::Point_2, typename Triangulation::Vertex_handle>::const_iterator found = vertices.find(edge.first);
      if (found == vertices.end()) {
        first_vertex = t.insert(edge.first, search_start);
        vertices[edge.first] = first_vertex;
      } else {
        first_vertex = found->second;
      } search_start = first_vertex->face();
      found = vertices.find(edge.second);
      if (found == vertices.end()) {
        second_vertex = t.insert(edge.second, search_start);
        vertices[edge.second] = second_vertex;
      } else {
        second_vertex = found->second;
      } search_start = second_vertex->face();
      edges_to_insert.emplace_back(first_vertex, second_vertex);
    }

    // Insert edges
    for (auto const& edge: edges_to_insert) {
      t.odd_even_insert_constraint(edge.first, edge.second);
    }
  }

  // Add edges of the polygon to the triangulation
  void add_to_triangulation_odd_even(const Multipolygon_with_holes_2<Kernel, PolygonContainer>& multipolygon) {
    std::unordered_set<std::pair<typename Kernel::Point_2, typename Kernel::Point_2>,
                       boost::hash<std::pair<typename Kernel::Point_2, typename Kernel::Point_2>>> unique_edges;

    // Get unique edges
    for (auto const& polygon: multipolygon.polygons()) {
      for (auto const& edge: polygon.outer_boundary().edges()) {
        if (edge.source() == edge.target()) continue;
        std::pair<typename Kernel::Point_2, typename Kernel::Point_2> pair = (edge.source() < edge.target())?
        std::make_pair(edge.source(), edge.target()) : std::make_pair(edge.target(), edge.source());
        auto inserted = unique_edges.insert(pair);
        if (!inserted.second) unique_edges.erase(inserted.first);
      }
      for (auto const& hole: polygon.holes()) {
        for (auto const& edge: hole.edges()) {
          if (edge.source() == edge.target()) continue;
          std::pair<typename Kernel::Point_2, typename Kernel::Point_2> pair = (edge.source() < edge.target())?
          std::make_pair(edge.source(), edge.target()) : std::make_pair(edge.target(), edge.source());
          auto inserted = unique_edges.insert(pair);
          if (!inserted.second) unique_edges.erase(inserted.first);
        }
      }
    }

    // Insert vertices
    std::unordered_map<typename Kernel::Point_2, typename Triangulation::Vertex_handle> vertices;
    std::vector<std::pair<typename Triangulation::Vertex_handle, typename Triangulation::Vertex_handle>> edges_to_insert;
    edges_to_insert.reserve(unique_edges.size());
    for (auto const& edge: unique_edges) {
      typename Triangulation::Vertex_handle first_vertex, second_vertex;
      typename std::unordered_map<typename Kernel::Point_2, typename Triangulation::Vertex_handle>::const_iterator found = vertices.find(edge.first);
      if (found == vertices.end()) {
        first_vertex = t.insert(edge.first, search_start);
        vertices[edge.first] = first_vertex;
      } else {
        first_vertex = found->second;
      } search_start = first_vertex->face();
      found = vertices.find(edge.second);
      if (found == vertices.end()) {
        second_vertex = t.insert(edge.second, search_start);
        vertices[edge.second] = second_vertex;
      } else {
        second_vertex = found->second;
      } search_start = second_vertex->face();
      edges_to_insert.emplace_back(first_vertex, second_vertex);
    }

    // Insert edges
    for (auto const& edge: edges_to_insert) {
      t.odd_even_insert_constraint(edge.first, edge.second);
    }
  }

  // Label a region of adjacent triangles without passing through constraints
  // adjacent triangles that involve passing through constraints are added to to_check
  void label_region(typename Triangulation::Face_handle face, int label,
                    std::list<typename Triangulation::Face_handle>& to_check,
                    std::list<int>& to_check_added_by) {
    // std::cout << "Labelling region with " << label << std::endl;
    std::list<typename Triangulation::Face_handle> to_check_in_region;
    face->label() = label;
    to_check_in_region.push_back(face);
    face->processed() = true; // processed means added to a list (to ensure elements are only added once)

    while (!to_check_in_region.empty()) {
      for (int neighbour = 0; neighbour < 3; ++neighbour) {
        if (!t.is_constrained(typename Triangulation::Edge(to_check_in_region.front(), neighbour))) {
          if (to_check_in_region.front()->neighbor(neighbour)->label() == 0) { // unlabelled
            to_check_in_region.front()->neighbor(neighbour)->label() = label;
            to_check_in_region.push_back(to_check_in_region.front()->neighbor(neighbour));
            to_check_in_region.front()->neighbor(neighbour)->processed() = true;
          }
        } else { // constrained
          if (!to_check_in_region.front()->neighbor(neighbour)->processed()) { // not added to to_check
            to_check.push_back(to_check_in_region.front()->neighbor(neighbour));
            to_check_added_by.push_back(label);
            to_check_in_region.front()->neighbor(neighbour)->processed() = true;
          }
        }
      } to_check_in_region.pop_front();
    }
  }

  // Label triangles in triangulation
  void label_triangulation_odd_even() {

    // Simplify collinear edges (gets rid of order dependency)
    for (auto vertex: t.all_vertex_handles()) {
      typename Triangulation::Edge_circulator first_edge = t.incident_edges(vertex);
      typename Triangulation::Edge_circulator current_edge = first_edge;
      std::vector<typename Triangulation::Edge> incident_constrained_edges;
      do {
        if (t.is_constrained(*current_edge)) {
          incident_constrained_edges.push_back(*current_edge);
        } ++current_edge;
      } while (current_edge != first_edge);
      if (incident_constrained_edges.size() == 2) {
        typename Kernel::Point_2 v1 = incident_constrained_edges.front().first->vertex(incident_constrained_edges.front().first->ccw(incident_constrained_edges.front().second))->point();
        typename Kernel::Point_2 v2 = incident_constrained_edges.back().first->vertex(incident_constrained_edges.back().first->ccw(incident_constrained_edges.back().second))->point();
        if (CGAL::collinear(v1, vertex->point(), v2)) {
          // std::cout << "Collinear points" << std::endl;
          // std::cout << "v1: " << v1 << std::endl;
          // std::cout << "in: " << vertex->point() << std::endl;
          // std::cout << "v2: " << v2 << std::endl;
          t.remove_incident_constraints(vertex);
          t.remove(vertex);
          t.insert_constraint(v1, v2);
        }
      }
    }

    // Init labels
    for (auto const face: t.all_face_handles()) {
      face->label() = 0;
      face->processed() = false;
    }

    // Label exterior with label -1, marking it as processed and
    // putting interior triangles adjacent to it in to_check
    std::list<typename Triangulation::Face_handle> to_check;
    std::list<int> to_check_added_by;
    label_region(t.infinite_face(), -1, to_check, to_check_added_by);

    // Label region of front element to_check list
    while (!to_check.empty()) {

      if (to_check.front()->label() == 0) { // label = 0 means not labelled yet
        if (to_check_added_by.front() < 0) {
          label_region(to_check.front(), number_of_polygons+1, to_check, to_check_added_by);
          ++number_of_polygons;
        } else {
          label_region(to_check.front(), -(number_of_holes+2), to_check, to_check_added_by);
          ++number_of_holes;
        }
      } to_check.pop_front();
      to_check_added_by.pop_front();

    } // std::cout << number_of_polygons << " polygons with " << number_of_holes << " holes in triangulation" << std::endl;
  }

  // Reconstruct ring boundary starting from an edge (face + opposite vertex) that is part of it
  void reconstruct_ring(std::list<typename Kernel::Point_2>& ring,
                        typename Triangulation::Face_handle face_adjacent_to_boundary,
                        int opposite_vertex) {
    // std::cout << "Reconstructing ring for face " << face_adjacent_to_boundary->label() << "..." << std::endl;

    // Create ring
    typename Triangulation::Face_handle current_face = face_adjacent_to_boundary;
    int current_opposite_vertex = opposite_vertex;
    do {
      CGAL_assertion(current_face->label() == face_adjacent_to_boundary->label());
      current_face->processed() = true;
      typename Triangulation::Vertex_handle pivot_vertex = current_face->vertex(current_face->cw(current_opposite_vertex));
      // std::cout << "\tAdding point " << pivot_vertex->point() << std::endl;
      ring.push_back(pivot_vertex->point());
      typename Triangulation::Face_circulator fc = t.incident_faces(pivot_vertex, current_face);
      do {
        ++fc;
      } while (fc->label() != current_face->label());
      current_face = fc;
      current_opposite_vertex = fc->cw(fc->index(pivot_vertex));
    } while (current_face != face_adjacent_to_boundary ||
             current_opposite_vertex != opposite_vertex);

    // Start at lexicographically smallest vertex
    typename std::list<typename Kernel::Point_2>::iterator smallest_vertex = ring.begin();
    for (typename std::list<typename Kernel::Point_2>::iterator current_vertex = ring.begin();
         current_vertex != ring.end(); ++current_vertex) {
      if (*current_vertex < *smallest_vertex) smallest_vertex = current_vertex;
    }
    if (ring.front() != *smallest_vertex) {
      ring.splice(ring.begin(), ring, smallest_vertex, ring.end());
    }
  }

  // Reconstruct multipolygon based on the triangles labeled as inside the polygon
  void reconstruct_multipolygon() {
    mp.clear();
    std::vector<Polygon_2<Kernel, PolygonContainer>> polygons; // outer boundaries
    std::vector<std::set<Polygon_2<Kernel, PolygonContainer>, Polygon_less>> holes; // holes are ordered (per polygon)
    polygons.resize(number_of_polygons);
    holes.resize(number_of_polygons);

    for (auto const face: t.all_face_handles()) {
      face->processed() = false;
    }
    for (auto const &face: t.finite_face_handles()) {
      if (face->label() < 1) continue; // exterior triangle
      if (face->processed()) continue; // already reconstructed
      for (int opposite_vertex = 0; opposite_vertex < 3; ++opposite_vertex) {
        if (face->label() == face->neighbor(opposite_vertex)->label()) continue; // not adjacent to boundary

        // Reconstruct ring
        std::list<typename Kernel::Point_2> ring;
        reconstruct_ring(ring, face, opposite_vertex);

        // Put ring in polygons
        Polygon_2<Kernel, PolygonContainer> polygon(ring.begin(), ring.end());
        // std::cout << "Reconstructed ring for polygon " << face->label() << " with ccw? " << (polygon.orientation() == CGAL::COUNTERCLOCKWISE) << std::endl;
        if (polygon.orientation() == CGAL::COUNTERCLOCKWISE) {
          polygons[face->label()-1] = polygon;
        } else {
          holes[face->label()-1].insert(polygon);
        } break;
      }
    }

    // Create polygons with holes and put in multipolygon
    std::set<Polygon_with_holes_2<Kernel, PolygonContainer>, Polygon_with_holes_less> ordered_polygons;
    for (int i = 0; i < polygons.size(); ++i) {
      ordered_polygons.insert(Polygon_with_holes_2<Kernel, PolygonContainer>(polygons[i], holes[i].begin(), holes[i].end()));
    }
    for (auto const& polygon: ordered_polygons) {
      // std::cout << "Adding polygon " << polygon << std::endl;
      mp.add_polygon(polygon);
    }
  }

  // Erases the triangulation.
  void clear() {
    t.clear();
  }

  /// @}

  /// \name Access Functions
  /// @{

  Triangulation& triangulation() {
    return t;
  }

  Multipolygon_with_holes_2<Kernel, PolygonContainer> multipolygon() {
    return mp;
  }

  /// @}


protected:
  Triangulation t;
  Multipolygon_with_holes_2<Kernel, PolygonContainer> mp;
  int number_of_polygons, number_of_holes;
  typename Triangulation::Face_handle search_start;
};

} // namespace Polygon_repair_2
} // namespace CGAL

#endif // CGAL_POLYGON_REPAIR_2_H