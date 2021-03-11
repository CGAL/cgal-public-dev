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

// Testing.
#include "../../../../../test/Levels_of_detail/include/Saver.h"
#include "../../../../../test/Levels_of_detail/include/Utilities.h"

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
    using Base_triangulation = typename Building::Base::Triangulation;
    using Triangulation = typename Base_triangulation::Delaunay;
    using Face_2 = internal::Partition_face_2<Traits>;
    using Face_3 = internal::Partition_edge_3<Traits>;
    using Polyhedron = internal::Partition_face_3<Traits>;
    using Vertex_handle = typename Triangulation::Vertex_handle;
    using Face_circulator = typename Triangulation::Face_circulator;
    using Face_handle = typename Triangulation::Face_handle;

    using Polygon = std::vector<Point_3>;
    using Indices = std::vector<std::size_t>;

    using Nearest_face_neighbor_query = internal::Nearest_face_neighbor_query<Traits>;
    using Coplanar_region = internal::Coplanar_region<Traits>;
    using Region_growing = internal::Region_growing<
    std::vector<Polygon>, Nearest_face_neighbor_query, Coplanar_region>;

    Building_builder(
      const Partition& partition,
      const FT distance_threshold) :
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

    void add_lod2_from_partition_3(
      Building& building) const {

      std::vector<Face_3> roofs;
      std::vector<Segment_3> segments;
      create_roofs(building.bottom_z, roofs);

      create_edges2(roofs, segments, building);
      create_base2(roofs, building);

      create_walls2(segments, building);
      create_roofs2(roofs, building);
    }

    void add_lod2_from_partition_2(
      const std::vector<Plane_3>& roof_planes,
      Building& building) const {

      create_edges2_from_edges1(building);
      create_base2_from_base1(building);

      Base_triangulation base;
      create_triangulation(roof_planes, base);

      create_roofs2_from_partition_2(base, building);
      create_walls2_from_partition_2(base, building);
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

    void create_edges2_from_edges1(
      Building& building) const {

      const auto& edges1 = building.edges1;
      CGAL_assertion(!edges1.empty());
      auto& edges2 = building.edges2;
      edges2 = edges1;
    }

    void create_edges2(
      const std::vector<Face_3>& roofs,
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

    void create_base2_from_base1(
      Building& building) const {

      const auto& base1 = building.base1;
      CGAL_assertion(!base1.empty());
      auto& base2 = building.base2;
      base2 = base1;
    }

    void create_base2(
      const std::vector<Face_3>& roofs,
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

    void create_walls2_from_partition_2(
      const Base_triangulation& base,
      Building& building) const {

      const FT bottom_z = building.bottom_z;
      auto& walls = building.walls2;
      walls.clear();

      std::vector<Polygon> polygons;

      const auto& tri = base.delaunay;
      for (auto fh = tri.finite_faces_begin();
      fh != tri.finite_faces_end(); ++fh) {

        for (std::size_t k = 0; k < 3; ++k) {
          const auto& fhn = fh->neighbor(k);

          const std::size_t i1 = (k + 1) % 3;
          const std::size_t j1 = (k + 2) % 3;

          const std::size_t nk = fhn->index(fh);
          const std::size_t i2 = (nk + 2) % 3;
          const std::size_t j2 = (nk + 1) % 3;

          create_wall2(bottom_z, i1, j1, i2, j2, fh, fhn, polygons);
        }
      }
      create_planar_items(polygons, true, walls);
    }

    void create_wall2(
      const FT bottom_z,
      const std::size_t i1, const std::size_t j1,
      const std::size_t i2, const std::size_t j2,
      const Face_handle& f1,
      const Face_handle& f2,
      std::vector<Polygon>& polygons) const {

      if (f1->info().label == f2->info().label)
        return;

      if (f1->info().label == std::size_t(-1))
        return;

      if (f2->info().label == std::size_t(-1)) {
        create_outer_wall2(bottom_z, i1, j1, f1, polygons);
        return;
      }

      if (f1->info().label != f2->info().label) {
        create_inner_wall2(bottom_z, i1, j1, i2, j2, f1, f2, polygons);
        return;
      }
    }

    void create_outer_wall2(
      const FT bottom_z,
      const std::size_t i,
      const std::size_t j,
      const Face_handle& f,
      std::vector<Polygon>& polygons) const {

      const auto& s = f->vertex(i)->point();
      const auto& t = f->vertex(j)->point();

      const FT sz = f->info().z[i];
      const FT tz = f->info().z[j];

      const Point_3 p1 = Point_3(s.x(), s.y(), bottom_z);
      const Point_3 p2 = Point_3(t.x(), t.y(), bottom_z);
      const Point_3 p3 = Point_3(t.x(), t.y(), tz);
      const Point_3 p4 = Point_3(s.x(), s.y(), sz);

      polygons.push_back({p1, p2, p3, p4});
    }

    void create_inner_wall2(
      const FT bottom_z,
      const std::size_t i1, const std::size_t j1,
      const std::size_t i2, const std::size_t j2,
      const Face_handle& f1,
      const Face_handle& f2,
      std::vector<Polygon>& polygons) const {

      const auto& s = f1->vertex(i1)->point();
      const auto& t = f1->vertex(j1)->point();

      const FT sz1 = f1->info().z[i1];
      const FT tz1 = f1->info().z[j1];

      const FT sz2 = f2->info().z[i2];
      const FT tz2 = f2->info().z[j2];

      if (
        CGAL::abs(sz1 - sz2) < m_distance_threshold &&
        CGAL::abs(tz1 - tz2) < m_distance_threshold)
        return;

      const Point_3 p1 = Point_3(s.x(), s.y(), bottom_z);
      const Point_3 p2 = Point_3(t.x(), t.y(), bottom_z);
      const Point_3 p3 = Point_3(t.x(), t.y(), tz1);
      const Point_3 p4 = Point_3(s.x(), s.y(), sz1);
      const Point_3 p5 = Point_3(t.x(), t.y(), tz2);
      const Point_3 p6 = Point_3(s.x(), s.y(), sz2);

      polygons.push_back({p1, p2, p3, p4});
      polygons.push_back({p1, p2, p5, p6});
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

    void create_roofs2_from_partition_2(
      const Base_triangulation& base,
      Building& building) const {

      auto& roofs = building.roofs2;
      roofs.clear();

      std::vector<Polygon> polygons;
      create_roof_polygons(base, polygons);
      create_planar_items(polygons, false, roofs);
    }

    void create_triangulation(
      const std::vector<Plane_3>& roof_planes,
      Base_triangulation& base) const {

      create_constrained_triangulation(base.delaunay);
      add_roof_labels(base.delaunay);
      update_heights(roof_planes, base.delaunay);

      save_triangulation(
        base, "/Users/monet/Documents/gf/lod/logs/buildings/tmp/delaunay");
    }

    void create_constrained_triangulation(
      Triangulation& tri) const {

      tri.clear();
      const auto& faces = m_partition.faces;
      for (const auto& face : faces) {

        const std::size_t label = face.label;
        const auto& segments = face.edges;
        const auto& neighbors = face.neighbors;

        for (std::size_t i = 0; i < segments.size(); ++i) {
          const auto& segment = segments[i];

          const auto& s = segment.source();
          const auto& t = segment.target();

          if (internal::are_equal_points_2(s, t))
            continue;

          const Vertex_handle svh = tri.insert(s);
          const Vertex_handle tvh = tri.insert(t);

          if (neighbors[i] >= 0) {
            const auto& face_neighbor = faces[neighbors[i]];
            const std::size_t label_neighbor = face_neighbor.label;

            if (label != label_neighbor && svh != tvh)
              tri.insert_constraint(svh, tvh);
          }
        }
      }
    }

    void add_roof_labels(
      Triangulation& tri) const {

      const auto& faces = m_partition.faces;
      for (auto fh = tri.finite_faces_begin();
      fh != tri.finite_faces_end(); ++fh) {

        const auto& p0 = fh->vertex(0)->point();
        const auto& p1 = fh->vertex(1)->point();
        const auto& p2 = fh->vertex(2)->point();

        const FT x = (p0.x() + p1.x() + p2.x()) / FT(3);
        const FT y = (p0.y() + p1.y() + p2.y()) / FT(3);

        const Point_2 b = Point_2(x, y);

        for (const auto& face : faces) {
          const auto& del = face.base.delaunay;

          if (del.number_of_vertices() < 3)
            continue;

          const auto bh = del.locate(b);
          if (!del.is_infinite(bh)) {
            fh->info().label = face.label;
            break;
          }
        }
      }
    }

    void update_heights(
      const std::vector<Plane_3>& roof_planes,
      Triangulation& tri) const {

      std::vector<Face_circulator> fhs;
      for (auto fh = tri.finite_faces_begin();
      fh != tri.finite_faces_end(); ++fh) {

        for (std::size_t k = 0; k < 3; ++k) {
          const auto& vh = fh->vertex(k);
          Face_circulator fc = tri.incident_faces(vh, fh);
          get_fhs(fc, fhs);
          const auto& p = vh->point();
          const FT z = update_height(p, roof_planes, fhs);
          fh->info().z[k] = z;
        }
      }
    }

    void get_fhs(
      Face_circulator& fc,
      std::vector<Face_circulator>& fhs) const {

      fhs.clear();
      if (fc.is_empty()) return;
      Face_circulator end = fc;
      do {
        fhs.push_back(fc); ++fc;
      } while (fc != end);
    }

    FT update_height(
      const Point_2& p,
      const std::vector<Plane_3>& roof_planes,
      const std::vector<Face_circulator>& fhs) const {

      std::map<std::size_t, std::size_t> pairs;
      for (const auto& fh : fhs) {
        const std::size_t label = fh->info().label;
        if (label != std::size_t(-1))
          pairs[label] = label;
      }

      if (pairs.size() == 0)
        return handle_0_neighbors(p, pairs, roof_planes);

      if (pairs.size() == 1)
        return handle_1_neighbor(p, pairs, roof_planes);

      if (pairs.size() == 2)
        return handle_2_neighbors(p, pairs, fhs[0], roof_planes);

      if (pairs.size() > 2)
        return handle_n_neighbors(p, pairs, fhs[0], roof_planes);

      return FT(0);
    }

    FT handle_0_neighbors(
      const Point_2& p,
      const std::map<std::size_t, std::size_t>& pairs,
      const std::vector<Plane_3>& roof_planes) const {

      return FT(0);
    }

    FT handle_1_neighbor(
      const Point_2& p,
      const std::map<std::size_t, std::size_t>& pairs,
      const std::vector<Plane_3>& roof_planes) const {

      const auto it = pairs.begin();
      const std::size_t label = it->second;
      const auto& plane = roof_planes[label];
      const Point_3 q = internal::position_on_plane_3(p, plane);

      return q.z();
    }

    FT handle_2_neighbors(
      const Point_2& p,
      const std::map<std::size_t, std::size_t>& pairs,
      const Face_handle& fh,
      const std::vector<Plane_3>& roof_planes) const {

      auto it = pairs.begin();
      const std::size_t label1 = it->second; ++it;
      const std::size_t label2 = it->second;

      const auto& plane1 = roof_planes[label1];
      const auto& plane2 = roof_planes[label2];

      const Point_3 q1 = internal::position_on_plane_3(p, plane1);
      const Point_3 q2 = internal::position_on_plane_3(p, plane2);

      const FT z1 = q1.z();
      const FT z2 = q2.z();

      /*
      if (CGAL::abs(z1 - z2) < m_distance_threshold)
        return CGAL::max(z1, z2); */

      if (CGAL::abs(z1 - z2) < m_distance_threshold)
        return (z1 + z2) / FT(2);

      const std::size_t label = fh->info().label;
      const auto& plane = roof_planes[label];
      const Point_3 q = internal::position_on_plane_3(p, plane);
      return q.z();
    }

    FT handle_n_neighbors(
      const Point_2& p,
      const std::map<std::size_t, std::size_t>& pairs,
      const Face_handle& fh,
      const std::vector<Plane_3>& roof_planes) const {

      std::vector< std::vector<std::size_t> > groups;
      create_groups(p, pairs, roof_planes, groups);

      std::vector<FT> zs;
      create_zs_avg(p, roof_planes, groups, zs);

      if (groups.size() == 0)
        return FT(0);

      if (groups.size() == 1)
        return zs[0];

      if (groups.size() > 1) {

        const std::size_t label = fh->info().label;
        for (std::size_t i = 0; i < groups.size(); ++i) {
          for (std::size_t j = 0; j < groups[i].size(); ++j) {
            if (label == groups[i][j])
              return zs[i];
          }
        }
      }

      return FT(0);
    }

    void create_zs_avg(
      const Point_2& p,
      const std::vector<Plane_3>& roof_planes,
      const std::vector< std::vector<std::size_t> >& groups,
      std::vector<FT>& zs) const {

      zs.clear();
      zs.reserve(groups.size());

      for (const auto& group : groups) {

        FT z = FT(0);
        for (const std::size_t label : group) {
          const auto& plane = roof_planes[label];
          const Point_3 q = internal::position_on_plane_3(p, plane);
          z += q.z();
        }
        z /= static_cast<FT>(group.size());
        zs.push_back(z);
      }
    }

    void create_zs_max(
      const Point_2& p,
      const std::vector<Plane_3>& roof_planes,
      const std::vector< std::vector<std::size_t> >& groups,
      std::vector<FT>& zs) const {

      zs.clear();
      zs.reserve(groups.size());

      for (const auto& group : groups) {

        FT z = -internal::max_value<FT>();
        for (const std::size_t label : group) {
          const auto& plane = roof_planes[label];
          const Point_3 q = internal::position_on_plane_3(p, plane);
          z = CGAL::max(z, q.z());
        }
        zs.push_back(z);
      }
    }

    void create_zs_min(
      const Point_2& p,
      const std::vector<Plane_3>& roof_planes,
      const std::vector< std::vector<std::size_t> >& groups,
      std::vector<FT>& zs) const {

      zs.clear();
      zs.reserve(groups.size());

      for (const auto& group : groups) {

        FT z = internal::max_value<FT>();
        for (const std::size_t label : group) {
          const auto& plane = roof_planes[label];
          const Point_3 q = internal::position_on_plane_3(p, plane);
          z = CGAL::min(z, q.z());
        }
        zs.push_back(z);
      }
    }

    void create_groups(
      const Point_2& p,
      const std::map<std::size_t, std::size_t>& pairs,
      const std::vector<Plane_3>& roof_planes,
      std::vector< std::vector<std::size_t> >& groups) const {

      groups.clear();

      std::vector<std::size_t> items;
      items.reserve(pairs.size());
      for (const auto& pair : pairs)
        items.push_back(pair.second);
      std::vector<bool> states(items.size(), false);

      std::vector<std::size_t> group;
      for (std::size_t i = 0; i < items.size(); ++i) {
        if (states[i]) continue;
        const std::size_t label1 = items[i];

        group.clear();
        group.push_back(label1);

        for (std::size_t j = 0; j < items.size(); ++j) {
          if (states[j]) continue;
          const std::size_t label2 = items[j];

          const auto& plane1 = roof_planes[label1];
          const auto& plane2 = roof_planes[label2];

          const Point_3 q1 = internal::position_on_plane_3(p, plane1);
          const Point_3 q2 = internal::position_on_plane_3(p, plane2);

          const FT z1 = q1.z();
          const FT z2 = q2.z();

          if (CGAL::abs(z1 - z2) < m_distance_threshold) {
            group.push_back(label2);
            states[j] = true;
          }
        }
        groups.push_back(group);
      }
    }

    void create_roof_polygons(
      const Base_triangulation& base,
      std::vector<Polygon>& polygons) const {

      const auto& tri = base.delaunay;
      for (auto fh = tri.finite_faces_begin();
      fh != tri.finite_faces_end(); ++fh) {
        if (fh->info().label == std::size_t(-1))
          continue;

        const auto& p0 = fh->vertex(0)->point();
        const auto& p1 = fh->vertex(1)->point();
        const auto& p2 = fh->vertex(2)->point();

        const FT z0 = fh->info().z[0];
        const FT z1 = fh->info().z[1];
        const FT z2 = fh->info().z[2];

        const Point_3 q0 = Point_3(p0.x(), p0.y(), z0);
        const Point_3 q1 = Point_3(p1.x(), p1.y(), z1);
        const Point_3 q2 = Point_3(p2.x(), p2.y(), z2);

        polygons.push_back({q0, q1, q2});
      }
    }

    void create_roofs2(
      const std::vector<Face_3>& faces,
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
      std::vector<Face_3>& roofs) const {

      // Create bounds.
      std::vector<Face_3> bounds;
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

    void create_bounds(std::vector<Face_3>& bounds) const {

      bounds.clear(); Face_3 face;
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
      const std::vector<Face_3>& faces) const {

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

    void save_triangulation(
      const Base_triangulation& base,
      const std::string path) const {

      const FT z = 0;
      std::size_t num_vertices = 0;
      internal::Indexer<Point_3> indexer;

      std::vector<Point_3> vertices;
      std::vector<Indices> faces;
      std::vector<CGAL::Color> fcolors;

      Polygon_inserter<Traits> inserter(faces, fcolors);
      auto output_vertices = std::back_inserter(vertices);
      auto output_faces = boost::make_function_output_iterator(inserter);

      base.output_with_label_color(
        indexer, num_vertices, output_vertices, output_faces, z);

      Saver<Traits> saver;
      saver.export_polygon_soup(vertices, faces, fcolors, path);
    }
  };

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_INTERNAL_BUILDING_BUILDER_H
