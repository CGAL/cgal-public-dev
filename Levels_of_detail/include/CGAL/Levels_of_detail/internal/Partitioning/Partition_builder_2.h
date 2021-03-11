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

#ifndef CGAL_LEVELS_OF_DETAIL_INTERNAL_PARTITION_BUILDER_2_H
#define CGAL_LEVELS_OF_DETAIL_INTERNAL_PARTITION_BUILDER_2_H

// STL includes.
#include <map>
#include <vector>

// CGAL includes.
#include <CGAL/property_map.h>

// Internal includes.
#include <CGAL/Levels_of_detail/internal/struct.h>

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

template<typename GeomTraits>
class Partition_builder_2 {

public:
  using Traits = GeomTraits;

  using FT = typename Traits::FT;
  using Point_2 = typename Traits::Point_2;
  using Segment_2 = typename Traits::Segment_2;

  using Triangulation = internal::Triangulation<Traits>;
  using Face_handle = typename Triangulation::Delaunay::Face_handle;
  using Vertex_handle = typename Triangulation::Delaunay::Vertex_handle;
  using Partition_edge_2 = internal::Partition_edge_2<Traits>;
  using Partition_face_2 = internal::Partition_face_2<Traits>;

  using Partition_2 = internal::Partition_2<Traits>;
  using LF_circulator = typename Triangulation::Delaunay::Line_face_circulator;

  Partition_builder_2(
    const bool with_visibility) :
  m_with_visibility(with_visibility)
  { }

  void build(
    const std::vector< std::vector<Segment_2> >& contours,
    Partition_2& partition_2) const {

    partition_2.clear();

    if (contours.size() == 0)
      return;

    // Create triangulation.
    Triangulation base;
    auto& tri = base.delaunay;
    tri.clear();

    std::vector<Vertex_handle> vhs;
    for (const auto& contour : contours) {

      vhs.clear();
      for (const auto& segment : contour)
        vhs.push_back(tri.insert(segment.source()));

      for (std::size_t i = 0; i < vhs.size(); ++i) {
        const std::size_t ip = (i + 1) % vhs.size();
        const auto vh1 = vhs[i];
        const auto vh2 = vhs[ip];

        if (vh1 != vh2)
          tri.insert_constraint(vh1, vh2);
      }
    }

    if (tri.number_of_faces() < 1)
      return;

    if (m_with_visibility)
      compute_visibility(base);

    std::map<Face_handle, int> fmap;
    create_faces(base, partition_2, fmap);
    create_face_neighbors(base, fmap, partition_2);
    create_edges(base, fmap, partition_2);
  }

  void build(
    const std::vector<Segment_2>& segments,
    Partition_2& partition_2) const {

    partition_2.clear();

    if (segments.size() < 3)
      return;

    // Create triangulation.
    Triangulation base;
    auto& tri = base.delaunay;
    tri.clear();

    for (const auto& segment : segments) {
      const auto svh = tri.insert(segment.source());
      const auto tvh = tri.insert(segment.target());
      if (svh != tvh)
        tri.insert_constraint(svh, tvh);
    }

    if (tri.number_of_faces() < 1)
      return;

    if (m_with_visibility)
      compute_visibility(base);

    std::map<Face_handle, int> fmap;
    create_faces(base, partition_2, fmap);
    create_face_neighbors(base, fmap, partition_2);
    create_edges(base, fmap, partition_2);
  }

private:
  const bool m_with_visibility;

  void create_faces(
    const Triangulation& base,
    Partition_2& partition_2,
    std::map<Face_handle, int>& fmap) const {

    const auto& tri = base.delaunay;

    partition_2.faces.clear();
    partition_2.faces.reserve(tri.number_of_faces());

    fmap.clear();

    Partition_face_2 pface;
    std::vector<Vertex_handle> vhs(3);

    int idx = 0;
    for (auto fit = tri.finite_faces_begin();
    fit != tri.finite_faces_end(); ++fit) {
      const Face_handle fh = static_cast<Face_handle>(fit);
      pface.base.delaunay.clear();

      for (std::size_t k = 0; k < 3; ++k)
        vhs[k] = pface.base.delaunay.insert(fh->vertex(k)->point());
      for (std::size_t k = 0; k < 3; ++k) {
        const std::size_t kp = (k + 1) % 3;
        if (vhs[k] != vhs[kp])
          pface.base.delaunay.insert_constraint(vhs[k], vhs[kp]);
      }

      if (fh->info().interior) {
        pface.inside = FT(1); pface.outside = FT(0);
        pface.visibility = Visibility_label::INSIDE;
      } else {
        pface.inside = FT(0); pface.outside = FT(1);
        pface.visibility = Visibility_label::OUTSIDE;
      }

      partition_2.faces.push_back(pface);
      fmap[fh] = idx; ++idx;
    }
  }

  void create_face_neighbors(
    const Triangulation& base,
    const std::map<Face_handle, int>& fmap,
    Partition_2& partition_2) const {

    const auto& tri = base.delaunay;

    int idx = 0;
    for (auto fit = tri.finite_faces_begin();
    fit != tri.finite_faces_end(); ++fit) {
      const Face_handle fh = static_cast<Face_handle>(fit);

      auto& edges = partition_2.faces[idx].edges;
      auto& neighbors = partition_2.faces[idx].neighbors;

      edges.clear(); edges.reserve(3);
      neighbors.clear(); neighbors.reserve(3);

      for (std::size_t k = 0; k < 3; ++k) {
        const Face_handle fhn = fh->neighbor(k);
        if (tri.is_infinite(fhn)) neighbors.push_back(-1);
        else {
          CGAL_assertion(fmap.find(fhn) != fmap.end());
          neighbors.push_back(fmap.at(fhn));
        }

        const auto& p1 = fh->vertex((k + 1) % 3)->point();
        const auto& p2 = fh->vertex((k + 2) % 3)->point();
        edges.push_back(Segment_2(p1, p2));
      }
      ++idx;
    }
  }

  void create_edges(
    const Triangulation& base,
    const std::map<Face_handle, int>& fmap,
    Partition_2& partition_2) const {

    const auto& tri = base.delaunay;

    partition_2.edges.clear();
    partition_2.edges.reserve(tri.number_of_faces());
    for (auto eh = tri.finite_edges_begin();
    eh != tri.finite_edges_end(); ++eh) {
      const Face_handle fh = eh->first;
      const std::size_t idx = eh->second;
      const Face_handle fhn = fh->neighbor(idx);

      const auto& p1 = fh->vertex((idx + 1) % 3)->point();
      const auto& p2 = fh->vertex((idx + 2) % 3)->point();

      int f1 = -1, f2 = -1;
      if (!tri.is_infinite(fh)) {
        if (fmap.find(fh) != fmap.end())
          f1 = fmap.at(fh);
      }
      if (!tri.is_infinite(fhn)) {
        if (fmap.find(fhn) != fmap.end())
          f2 = fmap.at(fhn);
      }
      partition_2.edges.push_back(Partition_edge_2(p1, p2, f1, f2));
    }
  }

  void compute_visibility(Triangulation& base) const {

    auto& tri = base.delaunay;

    // Bbox;
    std::vector<Point_2> points;
    for (auto fit = tri.finite_faces_begin();
    fit != tri.finite_faces_end(); ++fit) {
      const Face_handle fh = static_cast<Face_handle>(fit);
      const auto& p0 = fh->vertex(0)->point();
      const auto& p1 = fh->vertex(1)->point();
      const auto& p2 = fh->vertex(2)->point();
      points.push_back(p0);
      points.push_back(p1);
      points.push_back(p2);
    }

    std::vector<Point_2> bbox;
    CGAL::Identity_property_map<Point_2> pmap;
    internal::bounding_box_2(points, pmap, bbox);

    // Visibility.
    for (auto fit = tri.finite_faces_begin();
    fit != tri.finite_faces_end(); ++fit) {
      const Face_handle fh = static_cast<Face_handle>(fit);

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

      if (in > FT(1) / FT(2)) fh->info().interior = true;
      else fh->info().interior = false;
    }
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
};

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_INTERNAL_PARTITION_BUILDER_2_H
