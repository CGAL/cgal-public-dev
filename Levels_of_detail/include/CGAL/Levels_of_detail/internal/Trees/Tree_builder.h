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

#ifndef CGAL_LEVELS_OF_DETAIL_INTERNAL_TREE_BUILDER_H
#define CGAL_LEVELS_OF_DETAIL_INTERNAL_TREE_BUILDER_H

#include <CGAL/license/Levels_of_detail.h>

// STL includes.
#include <cmath>
#include <vector>

// Internal includes.
#include <CGAL/Levels_of_detail/internal/struct.h>

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

  template<typename GeomTraits>
  class Tree_builder {

  public:
    using Traits = GeomTraits;

    using FT = typename Traits::FT;
    using Point_2 = typename Traits::Point_2;
    using Vector_2 = typename Traits::Vector_2;
    using Segment_2 = typename Traits::Segment_2;

    using Tree_model = internal::Tree_model<Traits>;
    using Tree = internal::Tree<Traits>;
    using Edge = typename Tree::Edge;
    using Triangulation = typename Tree::Base::Triangulation::Delaunay;

    void add_lod0(
      const Tree_model& model, 
      const std::size_t min_faces_per_footprint,
      Tree& tree) const {

      std::vector<Edge>& edges = tree.edges0;
      Triangulation& tri = tree.base0.triangulation.delaunay;
      create_base01(
        model.center, model.radius, min_faces_per_footprint, 
        edges, tri);
    }

  private:
    void create_base01(
      const Point_2& center,
      const FT radius,
      const std::size_t min_faces_per_footprint,
      std::vector<Edge>& edges, Triangulation& tri) const {

      // Create boundary points.
      const std::size_t n = min_faces_per_footprint;
      std::vector<Point_2> points(n);
      for (std::size_t i = 0; i < n; ++i) {
        
        const FT angle = FT(2) * static_cast<FT>(CGAL_PI) * 
        (static_cast<FT>(i) / static_cast<FT>(n));

        points[i] = center + Vector_2(
          radius * static_cast<FT>(std::cos(CGAL::to_double(angle))),
          radius * static_cast<FT>(std::sin(CGAL::to_double(angle))));
      }

      // Create boundary edges.
      edges.clear(); edges.reserve(n); Edge edge;
      for (std::size_t i = 0; i < n; ++i) {
        const std::size_t ip = (i + 1) % n;
        const Point_2& s = points[i];
        const Point_2& t = points[ip];
        edge.segment = Segment_2(s, t);
        edges.push_back(edge);
      }
      CGAL_assertion(edges.size() == n);

      // Create base triangulation.
      tri.clear();
      for (const Edge& edge : edges) {
        const auto svh = tri.insert(edge.segment.source());
        const auto tvh = tri.insert(edge.segment.target());
        if (svh != tvh)
          tri.insert_constraint(svh, tvh);
      }
    }
  };

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_INTERNAL_TREE_BUILDER_H
