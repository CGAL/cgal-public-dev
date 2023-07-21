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
Multipolygon_with_holes_2<Kernel, PolygonContainer> repair(const Polygon_2<Kernel, PolygonContainer>& p) {
  CGAL::Polygon_repair_2::Polygon_repair_2<Kernel, PolygonContainer> pr;
  pr.add_to_triangulation(p);
  pr.label_triangulation();
  pr.reconstruct_multipolygon();
  return pr.multipolygon();
}

/// \ingroup PkgPolygonRepair2Functions
/// Repair a polygon with holes
template <class Kernel, class PolygonContainer>
Multipolygon_with_holes_2<Kernel, PolygonContainer> repair(const Polygon_with_holes_2<Kernel, PolygonContainer>& p) {
  CGAL::Polygon_repair_2::Polygon_repair_2<Kernel, PolygonContainer> pr;
  pr.add_to_triangulation(p);
  pr.label_triangulation();
  pr.reconstruct_multipolygon();
  return pr.multipolygon();
}

/// \ingroup PkgPolygonRepair2Functions
/// Repair a multipolygon with holes
template <class Kernel, class PolygonContainer>
Multipolygon_with_holes_2<Kernel, PolygonContainer> repair(const Multipolygon_with_holes_2<Kernel, PolygonContainer>& mp) {
  CGAL::Polygon_repair_2::Polygon_repair_2<Kernel, PolygonContainer> pr;
  pr.add_to_triangulation(mp);
  pr.label_triangulation();
  pr.reconstruct_multipolygon();
  return pr.multipolygon();
}

/*! \ingroup PkgPolygonRepair2Ref
 *
 * The class `Polygon_repair_2` builds on a constrained
 * triangulation to remove the parts of constraints that overlap an even number of times
 */
template <class Kernel = CGAL::Exact_predicates_inexact_constructions_kernel,
          class PolygonContainer = std::vector<typename Kernel::Point_2>>
class Polygon_repair_2 {
  int number_of_polygons, number_of_holes;
  std::vector<int> hole_nesting;
public:
  typedef CGAL::Triangulation_vertex_base_2<Kernel> Vertex_base;
  typedef CGAL::Constrained_triangulation_face_base_2<Kernel> Face_base;
  typedef CGAL::Triangulation_face_base_with_repair_info_2<Kernel, Face_base> Face_base_with_repair_info;
  typedef CGAL::Triangulation_data_structure_2<Vertex_base, Face_base_with_repair_info> Triangulation_data_structure;
  typedef CGAL::Exact_predicates_tag Tag; // assumed for now
  typedef CGAL::Constrained_Delaunay_triangulation_2<Kernel, Triangulation_data_structure, Tag> Constrained_Delaunay_triangulation;
  typedef Triangulation_with_odd_even_constraints_2<Constrained_Delaunay_triangulation> Triangulation;

  struct Polygon_compare {
    bool operator()(const Polygon_2<Kernel, PolygonContainer>& pa, const Polygon_2<Kernel, PolygonContainer>& pb) const {
      typename Polygon_2<Kernel, PolygonContainer>::Vertex_iterator va = pa.vertices_begin();
      typename Polygon_2<Kernel, PolygonContainer>::Vertex_iterator vb = pb.vertices_begin();
      while (va != pa.vertices_end() && vb != pb.vertices_end()) {
        if (*va != *vb) return *va < *vb;
        ++va;
        ++vb;
      } if (vb == pb.vertices_end()) return false;
      return true;
    }
  };

  struct Polygon_with_holes_compare {
    bool operator()(const Polygon_with_holes_2<Kernel, PolygonContainer>& pa, const Polygon_with_holes_2<Kernel, PolygonContainer>& pb) const {
      typename Polygon_2<Kernel, PolygonContainer>::Vertex_iterator va = pa.outer_boundary().vertices_begin();
      typename Polygon_2<Kernel, PolygonContainer>::Vertex_iterator vb = pb.outer_boundary().vertices_begin();
      while (va != pa.outer_boundary().vertices_end() && vb != pb.outer_boundary().vertices_end()) {
        if (*va != *vb) return *va < *vb;
        ++va;
        ++vb;
      } if (vb != pb.outer_boundary().vertices_end()) return true;
      else if (va != pa.outer_boundary().vertices_end()) return false;
      typename Polygon_with_holes_2<Kernel, PolygonContainer>::Hole_const_iterator ha = pa.holes_begin();
      typename Polygon_with_holes_2<Kernel, PolygonContainer>::Hole_const_iterator hb = pb.holes_begin();
      while (ha != pa.holes_end() && hb != pb.holes_end()) {
        va = ha->vertices_begin();
        vb = hb->vertices_begin();
        while (va != ha->vertices_end() && vb != hb->vertices_end()) {
          if (*va != *vb) return *va < *vb;
          ++va;
          ++vb;
        } if (vb != hb->vertices_end()) return true;
        else if (va != ha->vertices_end()) return false;
        ++ha;
        ++hb;
      } if (hb == pb.holes_end()) return false;
      return true;
    }
  };

  /// \name Creation
  Polygon_repair_2() : number_of_polygons(0), number_of_holes(0) {}

  /// \name Modifiers
  /// @{

  // Add edges of the polygon to the triangulation
  void add_to_triangulation(const Polygon_2<Kernel, PolygonContainer>& polygon) {
    for (auto const& edge: polygon.edges()) {
      t.odd_even_insert_constraint(edge.source(), edge.target());
    }
  }

  // Add edges of the polygon to the triangulation
  void add_to_triangulation(const Polygon_with_holes_2<Kernel, PolygonContainer>& polygon) {
    add_to_triangulation(polygon.outer_boundary());
    for (auto const& hole: polygon.holes()) {
       add_to_triangulation(hole);
    }
  }

  // Add edges of the polygon to the triangulation
  void add_to_triangulation(const Multipolygon_with_holes_2<Kernel, PolygonContainer>& multipolygon) {
    for (auto const& polygon: multipolygon.polygons()) {
      add_to_triangulation(polygon);
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
          if (to_check_in_region.front()->neighbor(neighbour)->label() == 0) {
            to_check_in_region.front()->neighbor(neighbour)->label() = label;
            to_check_in_region.push_back(to_check_in_region.front()->neighbor(neighbour));
            to_check_in_region.front()->neighbor(neighbour)->processed() = true;
          }
        } else {
          if (!to_check_in_region.front()->neighbor(neighbour)->processed()) {
            to_check.push_back(to_check_in_region.front()->neighbor(neighbour));
            to_check_added_by.push_back(label);
            to_check_in_region.front()->neighbor(neighbour)->processed() = true;
          }
        }
      } to_check_in_region.pop_front();
    }
  }

  // Label triangles in triangulation
  void label_triangulation() {
    for (auto const face: t.all_face_handles()) {
      face->label() = 0;
      face->processed() = false;
    }

    // Mark exterior as processed and put interior triangles adjacent to it in to_check
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
          hole_nesting.push_back(to_check_added_by.front()); // record nesting of current hole
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

    // Create ring
    typename Triangulation::Face_handle current_face = face_adjacent_to_boundary;
    int current_opposite_vertex = opposite_vertex;
    do {
      typename Triangulation::Vertex_handle pivot_vertex = current_face->vertex(current_face->cw(current_opposite_vertex));
      // std::cout << "Adding point " << pivot_vertex->point() << std::endl;
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
    } if (ring.front() != *smallest_vertex) {
      ring.splice(ring.begin(), ring, smallest_vertex, ring.end());
    }
  }

  // Reconstruct multipolygon based on the triangles labelled as inside the polygon
  void reconstruct_multipolygon() {
    mp.clear();
    std::vector<Polygon_2<Kernel, PolygonContainer>> polygons;
    std::vector<std::set<Polygon_2<Kernel, PolygonContainer>, Polygon_compare>> holes; // holes are ordered
    for (int i = 0; i < number_of_polygons; ++i) {
      polygons.emplace_back();
      holes.emplace_back();
    }

    for (auto const face: t.all_face_handles()) {
      face->processed() = false;
    } for (auto const &face: t.finite_face_handles()) {
      if (face->label() == -1) continue; // exterior triangle
      if (face->processed()) continue; // already reconstructed
      for (int opposite_vertex = 0; opposite_vertex < 3; ++opposite_vertex) {
        if (face->label() != face->neighbor(opposite_vertex)->label()) {

          std::list<typename Kernel::Point_2> ring;
          reconstruct_ring(ring, face, opposite_vertex);
          if (face->label() > 0) {
            polygons[face->label()-1].insert(polygons[face->label()-1].vertices_end(),
                                             ring.begin(), ring.end());
          } else {
            // std::cout << "Label: " << face->label() << " -> item " << -face->label()-2 << " -> in polygon " << hole_nesting[-face->label()-2] << std::endl;
            ring.push_back(ring.front());
            ring.pop_front();
            holes[hole_nesting[-face->label()-2]-1].insert(Polygon_2<Kernel, PolygonContainer>(ring.rbegin(), ring.rend()));
          }

          std::list<typename Triangulation::Face_handle> to_check;
          to_check.push_back(face);
          while (!to_check.empty()) {
            for (int neighbour = 0; neighbour < 3; ++neighbour) {
              if (to_check.front()->label() == to_check.front()->neighbor(neighbour)->label() &&
                  !to_check.front()->neighbor(neighbour)->processed()) {
                to_check.push_back(to_check.front()->neighbor(neighbour));
              }
            } to_check.front()->processed() = true;
            to_check.pop_front();
          }

          break;
        }
      }
    }

    // Create polygons with holes and put in multipolygon
    std::set<Polygon_with_holes_2<Kernel, PolygonContainer>, Polygon_with_holes_compare> ordered_polygons;
    for (int i = 0; i < polygons.size(); ++i) {
      ordered_polygons.insert(Polygon_with_holes_2<Kernel, PolygonContainer>(polygons[i], holes[i].begin(), holes[i].end()));
    } for (auto const& polygon: ordered_polygons) {
      std::cout << "Adding polygon " << polygon << std::endl;
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
};

} // namespace Polygon_repair_2
} // namespace CGAL

#endif // CGAL_POLYGON_REPAIR_2_H
