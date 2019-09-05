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
  using Segment_2 = typename Traits::Segment_2;

  using Triangulation = internal::Triangulation<Traits>;
  using Face_handle = typename Triangulation::Delaunay::Face_handle;
  using Vertex_handle = typename Triangulation::Delaunay::Vertex_handle;
  using Partition_edge_2 = internal::Partition_edge_2<Traits>;
  using Partition_face_2 = internal::Partition_face_2<Traits>;

  using Partition_2 = internal::Partition_2<Traits>;

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
      const auto vhs = tri.insert(segment.source());
      const auto vht = tri.insert(segment.target());
      if (vhs != vht)
        tri.insert_constraint(vhs, vht);
    }

    if (tri.number_of_faces() < 1)
      return;

    std::map<Face_handle, int> fmap;
    create_faces(base, partition_2, fmap);
    create_face_neighbors(base, fmap, partition_2);
    create_edges(base, fmap, partition_2);
  }

private:

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
};

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_INTERNAL_PARTITION_BUILDER_2_H
