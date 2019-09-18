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
#include <CGAL/Levels_of_detail/internal/utils.h>

// Spatial search.
#include <CGAL/Levels_of_detail/internal/Spatial_search/Nearest_face_neighbor_query.h>

// Shape detection.
#include <CGAL/Levels_of_detail/internal/Shape_detection/Region_growing.h>
#include <CGAL/Levels_of_detail/internal/Shape_detection/Coplanar_region.h>

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
    using Segment_3 = typename Traits::Segment_3;
    using Triangle_2 = typename Traits::Triangle_2;
    using Triangle_3 = typename Traits::Triangle_3;
    using Plane_3 = typename Traits::Plane_3;

    using Building = internal::Building<Traits>;
    using Boundary = internal::Boundary<Traits>;
    using Wall = internal::Building_wall<Traits>;
    using Roof = internal::Building_roof<Traits>;
    using Triangulation = typename Building::Base::Triangulation::Delaunay;
    using Face = internal::Partition_edge_3<Traits>;
    using Polyhedron = internal::Partition_face_3<Traits>;
    
    using Polygon = std::vector<Point_3>;
    using Indices = std::vector<std::size_t>;
    
    using Nearest_face_neighbor_query = internal::Nearest_face_neighbor_query<Traits>;
    using Coplanar_region = internal::Coplanar_region<Traits>;
    using Region_growing = internal::Region_growing<
    std::vector<Polygon>, Nearest_face_neighbor_query, Coplanar_region>;

    Building_builder(
      const Partition& partition,
      const FT distance_threshold = FT(0)) :
    m_partition(partition),
    m_angle_threshold(FT(10)),
    m_distance_threshold(distance_threshold)
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

    void add_lod2(
      Building& building) const {
    
      std::vector<Face> roofs;
      std::vector<Segment_3> segments;
      create_roofs(building.bottom_z, roofs);
      
      create_edges2(roofs, segments, building);
      create_base2(roofs, building);

      create_walls2(segments, building);
      create_roofs2(roofs, building);
    }

    void add_lod2_from_image(
      const std::vector<Plane_3>& roof_planes,
      Building& building) const {

      
    }

  private:
    const Partition& m_partition;
    const FT m_angle_threshold;
    const FT m_distance_threshold;

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
      const FT bottom_z = building.bottom_z;
      edges.clear();
      
      Boundary edge;
      for (const std::size_t fidx : findices) {
        
        const auto& face = m_partition.faces[fidx];
        const auto& neighbors = face.neighbors;
        
        const auto& nedges = face.edges;
        CGAL_assertion(nedges.size() == neighbors.size());
        
        // const auto& nconstr = face.constraints;
        // CGAL_assertion(nconstr.size() > 0);

        for (std::size_t i = 0; i < neighbors.size(); ++i) {
          const int nidx = neighbors[i];
          if (
            nidx < 0 || 
            // nconstr.at(nidx) || 
            m_partition.faces[nidx].visibility == Visibility_label::OUTSIDE) {
            
            edge.segment = nedges[i];
            edge.z = bottom_z; 
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

    void create_edges2(
      const std::vector<Face>& roofs,
      std::vector<Segment_3>& segments,
      Building& building) const {
      
      if (roofs.empty()) return;
      auto& edges = building.edges2;
      const FT bottom_z = building.bottom_z;

      Boundary edge;
      edges.clear(); segments.clear();
      
      for (std::size_t i = 0; i < roofs.size(); ++i) {
        const auto& polygon = roofs[i].polygon;
        for (std::size_t j = 0; j < polygon.size(); ++j) {
          const std::size_t jp = (j + 1) % polygon.size();

          const auto& p1 = polygon[j];
          const auto& p2 = polygon[jp];

          if (!is_interior_edge(p1, p2, i, roofs)) {
            segments.push_back(Segment_3(p1, p2));
            const Point_2 q1 = Point_2(p1.x(), p1.y());
            const Point_2 q2 = Point_2(p2.x(), p2.y());
            edge.segment = Segment_2(q1, q2);
            edge.z = bottom_z; 
            edges.push_back(edge);
          }
        }
      }
    }

    void create_base0(
      const std::vector<std::size_t>& findices,
      Building& building) const {

      auto& tri = building.base0.triangulation.delaunay;
      const std::size_t index = building.index;
      const FT z = building.bottom_z;

      // Create triangulation.
      tri.clear();
      for (const std::size_t fidx : findices) {
        const auto& face = m_partition.faces[fidx];
        const auto& nedges = face.edges;
        add_edges(nedges, z, tri);
      }

      // Update faces.
      update_faces(index, z, tri);      
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
          if (!base.is_infinite(handle)) {
            found = true;
            break;
          }
        }
        if (!found) {
          fh->info().interior = false;
          fh->info().tagged = false;
        }
      }
    }

    void create_base1(
      Building& building) const {
      
      const auto& base0 = building.base0;
      CGAL_assertion(!base0.empty());
      auto& base1 = building.base1;
      base1 = base0;
    }

    void create_base2(
      const std::vector<Face>& roofs,
      Building& building) const {
      
      if (roofs.empty()) return;

      const auto& edges = building.edges2;
      auto& tri = building.base2.triangulation.delaunay;
      const std::size_t index = building.index;
      const FT z = building.bottom_z;

      // Create triangulation.
      tri.clear();
      for (const auto& edge : edges) {
        const auto svh = tri.insert(edge.segment.source());
        const auto tvh = tri.insert(edge.segment.target());

        svh->info().z = z;
        tvh->info().z = z;

        if (svh != tvh)
          tri.insert_constraint(svh, tvh);
      }

      // Update faces.
      update_faces(index, z, tri);
      for (auto fh = tri.finite_faces_begin(); 
      fh != tri.finite_faces_end(); ++fh) {

        const Point_2 b = CGAL::barycenter(
          fh->vertex(0)->point(), FT(1),
          fh->vertex(1)->point(), FT(1),
          fh->vertex(2)->point(), FT(1));

        bool found = false;
        for (const auto& roof : roofs) {
          if (is_inside_polygon(b, roof.polygon)) {
            found = true;
            break;
          }
        }
        if (!found) {
          fh->info().interior = false;
          fh->info().tagged = false;
        }
      }
    }

    void create_walls1(
      Building& building) const {

      const auto& edges = building.edges1;
      const FT bottom_z = building.bottom_z;
      const FT top_z = building.top_z;
      auto& walls = building.walls1;
      
      std::vector<Polygon> polygons;
      polygons.resize(edges.size());
      for (std::size_t i = 0; i < edges.size(); ++i) {  
        const auto& edge = edges[i];
        const Point_2& s = edge.segment.source();
        const Point_2& t = edge.segment.target();
        
        const Point_3 p1 = Point_3(s.x(), s.y(), bottom_z);
        const Point_3 p2 = Point_3(t.x(), t.y(), bottom_z);
        const Point_3 p3 = Point_3(t.x(), t.y(), top_z);
        const Point_3 p4 = Point_3(s.x(), s.y(), top_z);
        
        polygons[i] = {p1, p2, p3, p4};
      }
      create_planar_items(polygons, true, walls);
    }

    void create_walls2(
      const std::vector<Segment_3>& segments,
      Building& building) const {

      if (segments.empty()) return;

      const FT bottom_z = building.bottom_z;
      auto& walls = building.walls2;

      std::vector<Polygon> polygons;
      polygons.resize(segments.size());
      for (std::size_t i = 0; i < segments.size(); ++i) {  
        const auto& segment = segments[i];
        const Point_3& s = segment.source();
        const Point_3& t = segment.target();
        
        const Point_3 p1 = Point_3(s.x(), s.y(), bottom_z);
        const Point_3 p2 = Point_3(t.x(), t.y(), bottom_z);
        const Point_3 p3 = Point_3(t.x(), t.y(), t.z());
        const Point_3 p4 = Point_3(s.x(), s.y(), s.z());
        
        polygons[i] = {p1, p2, p3, p4};
      }
      create_planar_items(polygons, true, walls);
    }

    void create_roofs1(
      Building& building) const {
    
      const auto& edges = building.edges1;
      const FT top_z = building.top_z;
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

        const Point_3 p1 = Point_3(a.x(), a.y(), top_z);
        const Point_3 p2 = Point_3(b.x(), b.y(), top_z);
        const Point_3 p3 = Point_3(c.x(), c.y(), top_z);

        roofs[0].triangles.push_back(Triangle_3(p1, p2, p3));
      }

      roofs[0].segments.reserve(edges.size());
      for (const auto& edge : edges) {
        const Point_2& p1 = edge.segment.source();
        const Point_2& p2 = edge.segment.target();

        const Point_3 a = Point_3(p1.x(), p1.y(), top_z);
        const Point_3 b = Point_3(p2.x(), p2.y(), top_z);
        
        roofs[0].segments.push_back(Segment_3(a, b));
      }
    }

    void create_roofs2(
      const std::vector<Face>& faces,
      Building& building) const {

      if (faces.empty()) return;
      auto& roofs = building.roofs2;

      std::vector<Polygon> polygons;
      polygons.reserve(faces.size());
      for (std::size_t i = 0; i < faces.size(); ++i) {
        const auto& polygon = faces[i].polygon;
        polygons.push_back(polygon);
      }
      create_planar_items(polygons, false, roofs);
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

    void update_faces(
      const std::size_t index,
      const FT z,
      Triangulation& tri) const {
      
      for (auto fh = tri.finite_faces_begin(); 
      fh != tri.finite_faces_end(); ++fh) {
        fh->info().urban_tag = Urban_object_type::BUILDING_ROOF;
        fh->info().object_index = index;
        fh->info().interior = true;
        fh->info().tagged = true;
        fh->info().z = {z, z, z};
      }
    }

    void create_roofs(
      const FT bottom_z,
      std::vector<Face>& roofs) const {

      // Create bounds.
      std::vector<Face> bounds; 
      create_bounds(bounds);

      // Find roofs.
      roofs.clear();
      for (const auto& face : bounds) {
        if (
        !internal::is_vertical_polygon(face.polygon, m_angle_threshold) &&
        !is_ground_polygon(bottom_z, face.polygon))
          roofs.push_back(face);
      }
    }

    bool is_ground_polygon(
      const FT bottom_z,
      const std::vector<Point_3>& polygon) const {

      for (const auto& p : polygon) {
        if (CGAL::abs(p.z() - bottom_z) > m_distance_threshold)
          return false;
      }
      return true;
    }

    void create_bounds(std::vector<Face>& bounds) const {

      bounds.clear(); Face face;
      const auto& polyhedrons = m_partition.faces;
      for (const auto& polyhedron : polyhedrons) {
        if (polyhedron.visibility == Visibility_label::OUTSIDE)
          continue;
        
        const auto& neighbors = polyhedron.neighbors;
        for (std::size_t i = 0; i < neighbors.size(); ++i) {
          const int idx = neighbors[i];
          if (idx >= 0 &&
              polyhedrons[idx].visibility == Visibility_label::OUTSIDE) {

            const auto& findices = polyhedron.faces[i];
            face.polygon.clear();
            for (const std::size_t fi : findices)
              face.polygon.push_back(polyhedron.vertices[fi]);
            bounds.push_back(face);
          }
        }
      }
    }

    bool is_inside_polygon(
      const Point_2& query, 
      const std::vector<Point_3>& poly_3) const {

      std::vector<Point_2> poly_2;
      internal::polygon_3_to_polygon_2(poly_3, poly_2);
      return internal::is_inside_polygon_2(query, poly_2);
    }

    bool is_interior_edge(
      const Point_3& p1, const Point_3& p2,
      const std::size_t curr_idx, 
      const std::vector<Face>& faces) const {

      for (std::size_t i = 0; i < faces.size(); ++i) {
        if (i == curr_idx) continue;
        
        const auto& polygon = faces[i].polygon;
        for (std::size_t j = 0; j < polygon.size(); ++j) {
          const std::size_t jp = (j + 1) % polygon.size();

          const auto& q1 = polygon[j];
          const auto& q2 = polygon[jp];

          if (internal::are_equal_edges_3(p1, p2, q1, q2))
            return true;
        }
      }
      return false;
    }

    template<typename Item>
    void create_planar_items(
      const std::vector<Polygon>& polygons,
      const bool vertical,
      std::vector<Item>& items) const {

      Nearest_face_neighbor_query neighbor_query(polygons);
      Coplanar_region region(polygons);
      Region_growing region_growing(
        polygons, neighbor_query, region);
      std::vector<Indices> regions;
      region_growing.detect(std::back_inserter(regions));

      items.clear(); items.reserve(regions.size());
      Indices neighbors; Item item;
      for (const auto& region : regions) {
        
        item.triangles.clear();
        item.segments.clear();
        for (const std::size_t idx : region) {
          const auto& polygon = polygons[idx];
          add_triangles(polygon, vertical, item.triangles);
          add_segments(polygons, idx, region, vertical, item.segments);
        }
        items.push_back(item);
      }
    }

    void add_triangles(
      const Polygon& polygon,
      const bool vertical,
      std::vector<Triangle_3>& triangles) const {

      if (vertical) { // do not remove this even if it is the same as below!
        triangles.push_back(
          Triangle_3(polygon[0], polygon[1], polygon[2]));
        triangles.push_back(
          Triangle_3(polygon[2], polygon[3], polygon[0]));
        return;
      }

      const auto& ref = polygon[0];
      for (std::size_t i = 1; i < polygon.size() - 1; ++i) {
        const std::size_t ip = i + 1;

        const auto& p1 = ref;
        const auto& p2 = polygon[i];
        const auto& p3 = polygon[ip];

        triangles.push_back(Triangle_3(p1, p2, p3));
      }
    }

    void add_segments(
      const std::vector<Polygon>& polygons,
      const std::size_t pidx,
      const Indices& indices,
      const bool vertical,
      std::vector<Segment_3>& segments) const {

      Indices skip;
      const auto& poly1 = polygons[pidx];
      for (std::size_t i = 0; i < poly1.size(); ++i) {
        const std::size_t ip = (i + 1) % poly1.size();  
        
        bool found = false;
        for (const std::size_t idx : indices) {
          if (idx == pidx) continue; 
          
          const auto& poly2 = polygons[idx];
          for (std::size_t j = 0; j < poly2.size(); ++j) {
            const std::size_t jp = (j + 1) % poly2.size();
            
            if (internal::are_equal_edges_3(
              poly1[i], poly1[ip], poly2[j], poly2[jp])) {
              skip.push_back(i);
              found = true; break;
            }
          }
          if (found) break;
        }
      }

      for (std::size_t i = 0; i < poly1.size(); ++i) {
        const std::size_t ip = (i + 1) % poly1.size();
        if (std::find(skip.begin(), skip.end(), i) != skip.end()) continue;
        if (vertical && internal::are_equal_points_2(
            Point_2(poly1[i].x(), poly1[i].y()),
            Point_2(poly1[ip].x(), poly1[ip].y())))
          segments.push_back(Segment_3(poly1[i], poly1[ip]));
        else if (!vertical)
          segments.push_back(Segment_3(poly1[i], poly1[ip]));
      }
    }
  };

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_INTERNAL_BUILDING_BUILDER_H
