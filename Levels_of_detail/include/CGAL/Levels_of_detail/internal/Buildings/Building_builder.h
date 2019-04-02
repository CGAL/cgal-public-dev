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

#ifndef CGAL_LEVELS_OF_DETAIL_INTERNAL_BUILDING_BUILDER_H
#define CGAL_LEVELS_OF_DETAIL_INTERNAL_BUILDING_BUILDER_H

#include <CGAL/license/Levels_of_detail.h>

// STL includes.
#include <cmath>
#include <vector>

// CGAL includes.
#include <CGAL/barycenter.h>
#include <CGAL/assertions.h>

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
  typename Partition_d,
  typename InputRange,
  typename PointMap>
  class Building_builder {

  public:
    using Traits = GeomTraits;
    using Partition = Partition_d;
    using Cluster = InputRange;
    using Point_map = PointMap;

    using FT = typename Traits::FT;
    using Point_2 = typename Traits::Point_2;
    using Point_3 = typename Traits::Point_3;
    using Segment_2 = typename Traits::Segment_2;
    using Triangle_2 = typename Traits::Triangle_2;
    using Triangle_3 = typename Traits::Triangle_3;

    using Building = internal::Building<Traits>;
    using Boundary = internal::Boundary<Traits>;
    using Wall = internal::Building_wall<Traits>;
    using Roof = internal::Building_roof<Traits>;
    using Triangulation = typename Building::Base::Triangulation::Delaunay;

    Building_builder(const Partition& partition) :
    m_partition(partition)
    { }

    void add_lod0(
      const std::vector<std::size_t>& base,
      const Cluster& cluster,
      const Point_map& point_map,
      Building& building) const {

      const bool success = create_bottom_z(cluster, point_map, building);
      if (!success) return;

      create_edges0(base, building);
      create_base0(base, building);
    }

    void add_lod1(
      const Extrusion_type extrusion_type,
      const Cluster& cluster,
      const Point_map& point_map,
      Building& building) const {

      const bool success = create_top_z(
        extrusion_type, cluster, point_map, building);
      if (!success) return;

      create_edges1(building);
      create_base1(building);

      create_walls1(building);
      create_roofs1(building);
    }

    void add_lod2() const {
    
      create_edges2();
      create_base2();

      create_walls2();
      create_roofs2();
    }

  private:
    const Partition& m_partition;

    bool create_bottom_z(
      const Cluster& cluster,
      const Point_map& point_map,
      Building& building) const {
    
      if (cluster.empty()) return false;
      const FT default_z = internal::max_value<FT>();

      building.bottom_z = default_z;
      for (const std::size_t idx : cluster) {
        const Point_3& p = get(point_map, idx);
        building.bottom_z = CGAL::min(building.bottom_z, p.z());
      }
      return true;
    }

    bool create_top_z(
      const Extrusion_type extrusion_type,
      const Cluster& cluster, 
      const Point_map& point_map, 
      Building& building) const {

      switch (extrusion_type) {
        case Extrusion_type::AVG: {
          return create_top_z_avg(cluster, point_map, building); }
        case Extrusion_type::MAX: {
          return create_top_z_max(cluster, point_map, building); }
        default: {
          return create_top_z_max(cluster, point_map, building); }
      }
    }

    bool create_top_z_avg(
      const Cluster& cluster, 
      const Point_map& point_map, 
      Building& building) const {

      if (cluster.empty()) return false;
      building.top_z = 0;
      for (const std::size_t idx : cluster) {
        const Point_3& p = get(point_map, idx);
        building.top_z += p.z();
      }
      building.top_z /= static_cast<FT>(cluster.size());
      return true;
    }

    bool create_top_z_max(
      const Cluster& cluster, 
      const Point_map& point_map, 
      Building& building) const {

      if (cluster.empty()) return false;
      const FT default_z = -internal::max_value<FT>();

      building.top_z = default_z;
      for (const std::size_t idx : cluster) {
        const Point_3& p = get(point_map, idx);
        building.top_z = CGAL::max(building.top_z, p.z());
      }
      return true;
    }

    void create_edges0(
      const std::vector<std::size_t>& findices,
      Building& building) const {

      auto& edges = building.edges0;
      edges.clear();
      
      Boundary edge;
      for (const std::size_t fidx : findices) {
        
        const auto& face = m_partition.faces[fidx];
        const auto& neighbors = face.neighbors;
        
        const auto& nedges = face.edges;
        CGAL_assertion(nedges.size() == neighbors.size());
        
        const auto& nconstr = face.constraints;
        CGAL_assertion(nconstr.size() > 0);

        for (std::size_t i = 0; i < neighbors.size(); ++i) {
          const int nidx = neighbors[i];
          if (
            nidx < 0 || 
            // nconstr.at(nidx) || 
            m_partition.faces[nidx].visibility == Visibility_label::OUTSIDE) {
            
            edge.segment = nedges[i];
            edge.z = building.bottom_z; 
            edges.push_back(edge);
          }
        }
      }
    }

    void create_edges1(
      Building& building) const {
      
      const auto& edges0 = building.edges0;
      CGAL_assertion(!edges0.empty());
      auto& edges1 = building.edges1;
      edges1 = edges0;
    }

    void create_edges2() const {
      
    }

    void create_base0(
      const std::vector<std::size_t>& findices,
      Building& building) const {

      create_base(
        findices, building.index, building.bottom_z, 
        building.base0.triangulation.delaunay);
    }

    void create_base1(
      Building& building) const {
      
      const auto& base0 = building.base0;
      CGAL_assertion(!base0.empty());
      auto& base1 = building.base1;
      base1 = base0;
    }

    void create_base2() const {
      
    }

    void create_walls1(
      Building& building) const {

      create_walls(
        building.edges1, building.bottom_z, building.top_z,
        building.walls1);
    }

    void create_walls2() const {

    }

    void create_roofs1(
      Building& building) const {
    
      auto& roofs = building.roofs1;
      CGAL_assertion(!building.base1.empty());
      const auto& tri = building.base1.triangulation.delaunay;

      roofs.clear();
      roofs.resize(1);
      roofs[0].triangles.reserve(tri.number_of_faces());

      for (auto fh = tri.finite_faces_begin();
      fh != tri.finite_faces_end(); ++fh) {
        if (!fh->info().tagged) continue;
        
        const Point_2& a = fh->vertex(0)->point();
        const Point_2& b = fh->vertex(1)->point();
        const Point_2& c = fh->vertex(2)->point();

        const Point_3 p1 = Point_3(a.x(), a.y(), building.top_z);
        const Point_3 p2 = Point_3(b.x(), b.y(), building.top_z);
        const Point_3 p3 = Point_3(c.x(), c.y(), building.top_z);

        const Triangle_3 tri = Triangle_3(p1, p2, p3);
        roofs[0].triangles.push_back(tri);
      }
    }

    void create_roofs2() const {

    }

    void create_base(
      const std::vector<std::size_t>& findices,
      const std::size_t index,
      const FT z,
      Triangulation& tri) const {
      
      // Create triangulation.
      tri.clear();
      for (const std::size_t fidx : findices) {
        const auto& face = m_partition.faces[fidx];
        const auto& nedges = face.edges;
        add_edges(nedges, z, tri);
      }

      // Update faces.
      for (auto fh = tri.finite_faces_begin(); 
      fh != tri.finite_faces_end(); ++fh) {
        fh->info().urban_tag = Urban_object_type::BUILDING;
        fh->info().object_index = index;
        fh->info().tagged = true;
        fh->info().z = {z, z, z};
      }

      for (auto fh = tri.finite_faces_begin(); 
      fh != tri.finite_faces_end(); ++fh) {

        const Point_2 b = CGAL::barycenter(
          fh->vertex(0)->point(), FT(1),
          fh->vertex(1)->point(), FT(1),
          fh->vertex(2)->point(), FT(1));

        bool found = false;
        for (const std::size_t fidx : findices) {
          const auto& face = m_partition.faces[fidx];
          const auto& base = face.base.delaunay;

          const auto handle = base.locate(b);
          if (!base.is_infinite(handle))
            found = true;
        }
        if (!found)
          fh->info().tagged = false;
      }
    }

    void add_edges(
      const std::vector<Segment_2>& edges,
      const FT z,
      Triangulation& tri) const {

      // Create base triangulation.
      for (const auto& edge : edges) {
        auto svh = tri.insert(edge.source());
        auto tvh = tri.insert(edge.target());

        svh->info().z = z;
        tvh->info().z = z;

        if (svh != tvh)
          tri.insert_constraint(svh, tvh);
      }
    }

    void create_walls(
      const std::vector<Boundary>& edges,
      const FT bottom_z,
      const FT top_z,
      std::vector<Wall>& walls) const {

      walls.clear();
      walls.resize(edges.size());

      for (std::size_t i = 0; i < edges.size(); ++i) {  
        const auto& edge = edges[i];
        
        const Point_2& s = edge.segment.source();
        const Point_2& t = edge.segment.target();
        
        const Point_3 p1 = Point_3(s.x(), s.y(), bottom_z);
        const Point_3 p2 = Point_3(t.x(), t.y(), bottom_z);
        const Point_3 p3 = Point_3(t.x(), t.y(), top_z);
        const Point_3 p4 = Point_3(s.x(), s.y(), top_z);
        
        const Triangle_3 tri1 = Triangle_3(p1, p2, p3);
        const Triangle_3 tri2 = Triangle_3(p3, p4, p1);

        walls[i].triangles.push_back(tri1);
        walls[i].triangles.push_back(tri2);
      }
    }
  };

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_INTERNAL_BUILDING_BUILDER_H
