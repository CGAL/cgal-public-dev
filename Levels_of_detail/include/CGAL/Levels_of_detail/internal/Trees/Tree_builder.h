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

// LOD includes.
#include <CGAL/Levels_of_detail/enum.h>

// Internal includes.
#include <CGAL/Levels_of_detail/internal/struct.h>
#include <CGAL/Levels_of_detail/internal/number_utils.h>

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

  template<
  typename GeomTraits,
  typename InputRange,
  typename PointMap>
  class Tree_builder {

  public:
    using Traits = GeomTraits;
    using Input_range = InputRange;
    using Point_map = PointMap;

    using FT = typename Traits::FT;
    using Point_2 = typename Traits::Point_2;
    using Point_3 = typename Traits::Point_3;
    using Vector_2 = typename Traits::Vector_2;
    using Segment_2 = typename Traits::Segment_2;
    using Triangle_3 = typename Traits::Triangle_3;

    using Tree_model = internal::Tree_model<Traits>;
    using Tree = internal::Tree<Traits>;
    using Edge = typename Tree::Edge;
    using Iterator = typename Input_range::const_iterator;

    void add_lod0(
      const Tree_model& model,
      const std::vector<Iterator>& cluster,
      const Point_map& point_map,
      const std::size_t min_faces_per_footprint,  
      Tree& tree) const {

      create_bottom_z(cluster, point_map, tree);
      const std::size_t n = min_faces_per_footprint;
      create_edges0(model, n, tree);
      create_base0(model, tree);
    }

    void add_lod1(
      const Extrusion_type extrusion_type,
      const std::vector<Iterator>& cluster,
      const Point_map& point_map,
      Tree& tree) const {

      create_top_z(extrusion_type, cluster, point_map, tree);
      create_edges1(tree);
      create_base1(tree);
      create_trunk1(tree);
      create_crown1(tree);
    }

    void add_lod2() const {
    
    }

  private:
    void create_bottom_z(
      const std::vector<Iterator>& cluster, 
      const Point_map& point_map, 
      Tree& tree) const {

      CGAL_assertion(!cluster.empty());
      const FT default_z = internal::max_value<FT>();

      tree.bottom_z = default_z;
      for (const auto& it : cluster) {
        const Point_3& p = get(point_map, *it);
        tree.bottom_z = CGAL::min(tree.bottom_z, p.z());
      }
    }

    void create_top_z(
      const Extrusion_type extrusion_type,
      const std::vector<Iterator>& cluster, 
      const Point_map& point_map, 
      Tree& tree) const {

      switch (extrusion_type) {
        case Extrusion_type::AVG: {
          create_top_z_avg(cluster, point_map, tree);
          return; }
        case Extrusion_type::MAX: {
          create_top_z_max(cluster, point_map, tree);
          return; }
        default: {
          create_top_z_max(cluster, point_map, tree);
          return; }
      }
    }

    void create_top_z_avg(
      const std::vector<Iterator>& cluster, 
      const Point_map& point_map, 
      Tree& tree) const {

      CGAL_assertion(!cluster.empty());
      tree.top_z = 0;
      for (const auto& it : cluster) {
        const Point_3& p = get(point_map, *it);
        tree.top_z += p.z();
      }
      tree.top_z /= static_cast<FT>(cluster.size());
    }

    void create_top_z_max(
      const std::vector<Iterator>& cluster, 
      const Point_map& point_map, 
      Tree& tree) const {

      CGAL_assertion(!cluster.empty());
      const FT default_z = -internal::max_value<FT>();

      tree.top_z = default_z;
      for (const auto& it : cluster) {
        const Point_3& p = get(point_map, *it);
        tree.top_z = CGAL::max(tree.top_z, p.z());
      }
    }

    void create_edges0(
      const Tree_model& model, 
      const std::size_t n, 
      Tree& tree) const {
      
      auto& edges = tree.edges0;

      // Create boundary points.
      std::vector<Point_2> points(n);
      for (std::size_t i = 0; i < n; ++i) {
        const FT angle = FT(2) * static_cast<FT>(CGAL_PI) * 
        (static_cast<FT>(i) / static_cast<FT>(n));

        points[i] = model.center + Vector_2(
          model.radius * static_cast<FT>(std::cos(CGAL::to_double(angle))),
          model.radius * static_cast<FT>(std::sin(CGAL::to_double(angle))));
      }

      // Create boundary edges.
      edges.clear(); edges.reserve(n); Edge edge;
      for (std::size_t i = 0; i < n; ++i) {
        const std::size_t ip = (i + 1) % n;
        const Point_2& s = points[i];
        const Point_2& t = points[ip];
        edge.segment = Segment_2(s, t);
        edge.z = tree.bottom_z;
        edges.push_back(edge);
      }
      CGAL_assertion(edges.size() == n);
    }

    void create_edges1(Tree& tree) const {
      const auto& edges0 = tree.edges0;
      CGAL_assertion(!edges0.empty());
      auto& edges1 = tree.edges1;
      edges1 = edges0;
    }

    void create_base0(
      const Tree_model& model,  
      Tree& tree) const {
      
      const auto& edges = tree.edges0;
      auto& tri = tree.base0.triangulation.delaunay;

      // Create base triangulation.
      tri.clear();
      for (const Edge& edge : edges) {
        auto svh = tri.insert(edge.segment.source());
        auto tvh = tri.insert(edge.segment.target());

        svh->info().z = tree.bottom_z;
        tvh->info().z = tree.bottom_z;

        if (svh != tvh)
          tri.insert_constraint(svh, tvh);
      }

      // Update faces.
      for (auto fh = tri.finite_faces_begin(); 
      fh != tri.finite_faces_end(); ++fh) {
        fh->info().urban_tag = Urban_object_type::TREE;
        fh->info().object_index = tree.index;
        fh->info().tagged = true;
        fh->info().z = {tree.bottom_z, tree.bottom_z, tree.bottom_z};
      }
    }

    void create_base1(Tree& tree) const {
      const auto& base0 = tree.base0;
      CGAL_assertion(!base0.empty());
      auto& base1 = tree.base1;
      base1 = base0;
    }

    void create_trunk1(Tree& tree) const {
    
      auto& trunk = tree.trunk1;
      std::vector<Triangle_3>& triangles = trunk.triangles;
      CGAL_assertion(!tree.edges1.empty());
      const auto& edges = tree.edges1;

      triangles.clear();
      triangles.reserve(edges.size() * 2);

      for (const auto& edge : edges) {  
        const Segment_2& seg = edge.segment;

        const Point_2& s = seg.source();
        const Point_2& t = seg.target();
        
        const Point_3 p1 = Point_3(s.x(), s.y(), tree.bottom_z);
        const Point_3 p2 = Point_3(t.x(), t.y(), tree.bottom_z);
        const Point_3 p3 = Point_3(t.x(), t.y(), tree.top_z);
        const Point_3 p4 = Point_3(s.x(), s.y(), tree.top_z);
        
        const Triangle_3 tri1 = Triangle_3(p1, p2, p3);
        const Triangle_3 tri2 = Triangle_3(p3, p4, p1);

        triangles.push_back(tri1);
        triangles.push_back(tri2);
      }
    }

    void create_crown1(Tree& tree) const {
    
      auto& crown = tree.crown1;
      std::vector<Triangle_3>& triangles = crown.triangles;
      CGAL_assertion(!tree.base1.empty());
      const auto& tri = tree.base1.triangulation.delaunay;

      triangles.clear();
      triangles.reserve(tri.number_of_faces());

      for (auto fh = tri.finite_faces_begin();
      fh != tri.finite_faces_end(); ++fh) {
        
        const Point_2& a = fh->vertex(0)->point();
        const Point_2& b = fh->vertex(1)->point();
        const Point_2& c = fh->vertex(2)->point();

        const Point_3 p1 = Point_3(a.x(), a.y(), tree.top_z);
        const Point_3 p2 = Point_3(b.x(), b.y(), tree.top_z);
        const Point_3 p3 = Point_3(c.x(), c.y(), tree.top_z);

        const Triangle_3 tri = Triangle_3(p1, p2, p3);
        triangles.push_back(tri);
      }
    }
  };

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_INTERNAL_TREE_BUILDER_H
