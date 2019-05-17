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

  using Segment_2 = typename Traits::Segment_2;

  using Triangulation = internal::Triangulation<Traits>;
  using Face_handle = typename Triangulation::Delaunay::Face_handle;
  using Vertex_handle = typename Triangulation::Delaunay::Vertex_handle;
  using Partition_edge_2 = internal::Partition_edge_2<Traits>;
  using Partition_face_2 = internal::Partition_face_2<Traits>;

  using Partition_2 = internal::Partition_2<Traits>;

  Partition_builder_2(
    const std::vector<Segment_2>& segments) :
  m_segments(segments)
  { }

  void build(Partition_2& partition_2) const {

    // Create triangulation.
    Triangulation base;
    auto& tri = base.delaunay;
    tri.clear();

    for (const auto& segment : m_segments) {
      const auto vhs = tri.insert(segment.source());
      const auto vht = tri.insert(segment.target());
      if (vhs != vht)
        tri.insert_constraint(vhs, vht);
    }
    partition_2.clear();

    // Create faces.
    partition_2.faces.clear();
    partition_2.faces.reserve(tri.number_of_faces());
    std::map<Face_handle, int> fmap;

    Partition_face_2 pface;
    std::vector<Vertex_handle> vhs(3);

    int idx = 0;
    for (auto fh = tri.finite_faces_begin();
    fh != tri.finite_faces_end(); ++fh) {
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

    // Create face neighbors and edges.
    idx = 0;
    for (auto fh = tri.finite_faces_begin();
    fh != tri.finite_faces_end(); ++fh) {
      
      auto& edges = partition_2.faces[idx].edges; 
      auto& neighbors = partition_2.faces[idx].neighbors;

      edges.clear(); edges.reserve(3);
      neighbors.clear(); neighbors.reserve(3);

      for (std::size_t k = 0; k < 3; ++k) {
        const auto fhn = fh->neighbor(k);
        if (tri.is_infinite(fhn)) neighbors.push_back(-1);
        else neighbors.push_back(fmap.at(fhn));

        const auto& p1 = fh->vertex((k + 1) % 3)->point();
        const auto& p2 = fh->vertex((k + 2) % 3)->point();
        edges.push_back(Segment_2(p1, p2));
      }
      ++idx;
    }

    // Create edges.
    partition_2.edges.clear();
    partition_2.edges.reserve(tri.number_of_faces());
    for (auto eh = tri.finite_edges_begin(); 
    eh != tri.finite_edges_end(); ++eh) {
      const auto fh = eh->first;
      const std::size_t idx = eh->second;
      const auto fhn = fh->neighbor(idx);

      const auto& p1 = fh->vertex((idx + 1) % 3)->point();
      const auto& p2 = fh->vertex((idx + 2) % 3)->point();

      int f1 = -1, f2 = -1;
      if (!tri.is_infinite(fh)) f1 = fmap.at(fh);
      if (!tri.is_infinite(fhn)) f2 = fmap.at(fhn);
      partition_2.edges.push_back(Partition_edge_2(p1, p2, f1, f2));
    }
  }
private:
  const std::vector<Segment_2>& m_segments;
};

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_INTERNAL_PARTITION_BUILDER_2_H
