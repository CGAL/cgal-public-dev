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

#ifndef CGAL_LEVELS_OF_DETAIL_INTERNAL_WALLS_ESTIMATOR_H
#define CGAL_LEVELS_OF_DETAIL_INTERNAL_WALLS_ESTIMATOR_H

// STL includes.
#include <vector>
#include <utility>

// CGAL includes.
#include <CGAL/assertions.h>

// Internal includes.
#include <CGAL/Levels_of_detail/internal/utils.h>
#include <CGAL/Levels_of_detail/internal/struct.h>

// Spatial search.
#include <CGAL/Levels_of_detail/internal/Spatial_search/Nearest_face_neighbor_query.h>

// Shape detection.
#include <CGAL/Levels_of_detail/internal/Shape_detection/Region_growing.h>
#include <CGAL/Levels_of_detail/internal/Shape_detection/Coplanar_region.h>

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

  template<typename GeomTraits>
  class Building_walls_estimator {

  public:
    using Traits = GeomTraits;

    using FT = typename Traits::FT;
    using Point_2 = typename Traits::Point_2;
    using Point_3 = typename Traits::Point_3;
    using Vector_3 = typename Traits::Vector_3;
    using Plane_3 = typename Traits::Plane_3;

    using Indices = std::vector<std::size_t>;
    using Boundary = internal::Boundary<Traits>;
    using Approximate_face = internal::Partition_edge_3<Traits>;
    using Polygon = std::vector<Point_3>;

    using Nearest_face_neighbor_query = internal::Nearest_face_neighbor_query<Traits>;
    using Coplanar_region = internal::Coplanar_region<Traits>;
    using Region_growing = internal::Region_growing<
    std::vector<Polygon>, Nearest_face_neighbor_query, Coplanar_region>;
    using Triangulation = internal::Triangulation<Traits>;
    using CDT = typename Triangulation::Delaunay;
    using Vertex_handle = typename CDT::Vertex_handle;
    using Face_handle = typename CDT::Face_handle;
    using Edge = typename CDT::Edge;

    using Vhs = std::vector<Vertex_handle>;
    using Vh_pair = std::pair<Vertex_handle, Vertex_handle>;

    Building_walls_estimator(
      const std::vector<Boundary>& boundaries,
      const FT bottom_z,
      const FT top_z) :
    m_boundaries(boundaries),
    m_bottom_z(bottom_z),
    m_top_z(top_z),
    m_max_num_iters(100)
    { }

    void estimate(
      std::vector<Approximate_face>& walls) const {

      if (m_boundaries.empty())
        return;

      std::vector<Polygon> faces;
      faces.clear(); faces.reserve(m_boundaries.size());

      Polygon face;
      for (const auto& boundary : m_boundaries) {
        estimate_wall(boundary, face);
        faces.push_back(face);
      }
      CGAL_assertion(faces.size() == m_boundaries.size());

      std::vector<Indices> regions;
      detect_coplanar_walls(faces, regions);
      merge_coplanar_walls(faces, regions, walls);
    }

    void estimate_wall(
      const Boundary& boundary,
      Polygon& face) const {

      const Point_2& s = boundary.segment.source();
      const Point_2& t = boundary.segment.target();

      face.clear();
      face.reserve(4);

      face.push_back(Point_3(s.x(), s.y(), m_bottom_z));
      face.push_back(Point_3(t.x(), t.y(), m_bottom_z));
      face.push_back(Point_3(t.x(), t.y(), m_top_z));
      face.push_back(Point_3(s.x(), s.y(), m_top_z));
    }

  private:
    const std::vector<Boundary>& m_boundaries;
    const FT m_bottom_z;
    const FT m_top_z;
    const std::size_t m_max_num_iters;

    void detect_coplanar_walls(
      const std::vector<Polygon>& faces,
      std::vector<Indices>& regions) const {

      Nearest_face_neighbor_query neighbor_query(faces);
      Coplanar_region region(faces);

      Region_growing region_growing(
        faces,
        neighbor_query,
        region);

      regions.clear();
      region_growing.detect(std::back_inserter(regions));
    }

    void merge_coplanar_walls(
      const std::vector<Polygon>& faces,
      const std::vector<Indices>& regions,
      std::vector<Approximate_face>& walls) const {

      walls.clear();
      if (regions.empty()) return;
      walls.reserve(regions.size());

      for (const auto& region : regions) {
        const bool success = merge_faces(faces, region, walls);
        if (!success) {

          // Safety feature.
          for (const std::size_t idx : region) {
            const auto& polygon = faces[idx];
            const auto& p1 = polygon[0];
            const auto& p2 = polygon[1];

            if (internal::distance(p1, p2) > internal::tolerance<FT>()) {
              Approximate_face wall;
              wall.polygon = polygon;
              walls.push_back(wall);
            }
          }
        }
      }
    }

    bool merge_faces(
      const std::vector<Polygon>& faces,
      const Indices& region,
      std::vector<Approximate_face>& merged) const {

      CGAL_assertion(region.size() > 0);

      Vector_3 n3; Point_3 b3;
      bool success = internal::compute_normal_3(faces[region[0]], n3);
      if (!success) return false;
      internal::compute_barycenter_3(faces, region, b3);
      const Plane_3 plane = Plane_3(b3, n3);

      std::vector<Polygon> rotated;
      rotated.reserve(region.size());
      Polygon poly;
      for (const std::size_t idx : region) {
        poly.clear();
        const auto& points_3 = faces[idx];
        for (const auto& p3 : points_3) {
          const auto& p2 = internal::to_2d(p3, b3, plane);
          poly.push_back(Point_3(p2.x(), p2.y(), FT(0)));
        }
        rotated.push_back(poly);
      }

      CDT cdt;
      triangulate(rotated, cdt);
      if (cdt.number_of_faces() == 0)
        return false;

      Polygon merged_face;
      success = create_merged_face(cdt, merged_face);
      if (!success || merged_face.size() < 3) return false;
      fix_orientation(merged_face);

      Approximate_face face;
      face.polygon.reserve(merged_face.size());
      for (const auto& p2 : merged_face) {
        const auto& q = Point_2(p2.x(), p2.y());
        face.polygon.push_back(internal::to_3d(q, b3, plane));
      }
      merged.push_back(face);
      return true;
    }

    void triangulate(
      const std::vector<Polygon>& faces,
      CDT& cdt) const {

			std::vector<Vhs> vhs;
      std::vector<Vh_pair> updated_vhs;
      insert_points(faces, cdt, vhs);
      update_constraints(faces, vhs, updated_vhs);
      insert_constraints(updated_vhs, cdt);
    }

    void insert_points(
      const std::vector<Polygon>& faces,
      CDT& cdt,
      std::vector<Vhs>& vhs) const {

      CGAL_assertion(faces.size() > 0);
      cdt.clear(); vhs.clear();
      vhs.resize(faces.size());
			for (std::size_t i = 0; i < faces.size(); ++i) {
				const auto& face = faces[i];

				vhs[i].resize(face.size());
				for (std::size_t j = 0; j < face.size(); ++j) {
          const auto& p = face[j];
					vhs[i][j] = cdt.insert(Point_2(p.x(), p.y()));
					vhs[i][j]->info().z = p.z();
				}
			}
    }

    void update_constraints(
      const std::vector<Polygon>& faces,
      const std::vector<Vhs>& vhs,
      std::vector<Vh_pair>& updated_vhs) const {

      CGAL_assertion(faces.size() > 0);
      for (std::size_t i = 0; i < faces.size(); ++i) {
        for (std::size_t j = 0; j < faces[i].size(); ++j) {
          const std::size_t jp = (j + 1) % faces[i].size();
          if (is_boundary_edge(
            faces[i][j], faces[i][jp], i, faces)) {

            const auto updated_vh = std::make_pair(vhs[i][j], vhs[i][jp]);
            updated_vhs.push_back(updated_vh);
          }
        }
      }
    }

    bool is_boundary_edge(
      const Point_3& p1, const Point_3& p2,
      const std::size_t face_index,
      const std::vector<Polygon>& faces) const {

      CGAL_assertion(faces.size() > 0);
      for (std::size_t i = 0; i < faces.size(); ++i) {
        if (i == face_index) continue;
        for (std::size_t j = 0; j < faces[i].size(); ++j) {
          const std::size_t jp = (j + 1) % faces[i].size();
          if (internal::are_equal_edges_3(
            p1, p2, faces[i][j], faces[i][jp]))
            return false;
        }
      }
      return true;
    }

    void insert_constraints(
      const std::vector<Vh_pair>& updated_vhs,
      CDT& cdt) const {

      CGAL_assertion(updated_vhs.size() > 0);
      for (const auto& vh : updated_vhs) {
        if (vh.first != vh.second)
          cdt.insert_constraint(vh.first, vh.second);
      }
    }

    bool create_merged_face(
      const CDT& cdt,
      Polygon& merged_face) const {

      merged_face.clear();
      if (cdt.number_of_faces() == 0) return false;
      Face_handle fh;
      bool success = find_first_face_handle(cdt, fh);
      if (!success) return false;
      success = traverse_cdt(fh, cdt, merged_face);
      if (!success) return false;
      return true;
    }

    bool find_first_face_handle(
      const CDT& cdt,
      Face_handle& fh) const {

      for (auto fit = cdt.finite_faces_begin();
      fit != cdt.finite_faces_end(); ++fit) {
        fh = static_cast<Face_handle>(fit);

        const auto& vh1 = fh->vertex(0);
        const auto& vh2 = fh->vertex(1);
        const auto& vh3 = fh->vertex(2);

        const auto& p1 = vh1->point();
        const auto& p2 = vh2->point();
        const auto& p3 = vh3->point();

        for (std::size_t k = 0; k < 3; ++k) {
          const Edge edge = std::make_pair(fh, k);
          if (cdt.is_constrained(edge))
            return true;
        }
      }
      return false;
    }

    bool traverse_cdt(
      const Face_handle& fh,
      const CDT& cdt,
      Polygon& face) const {

      Edge edge; face.clear();
      const bool success = find_first_edge(fh, cdt, edge);
      if (!success) return false;

      CGAL_assertion(edge.second >= 0 && edge.second <= 2);
      auto vh = edge.first->vertex((edge.second + 2) % 3);
      auto end = vh;
      if (vh->info().z == vh->info().default_z)
        return false;

      const auto& p = vh->point();
      face.push_back(Point_3(p.x(), p.y(), vh->info().z));
      std::size_t num_iters = 0;
      do {
        get_next_vertex_handle(cdt, vh, edge);
        const auto& q = vh->point();
        if (vh->info().z == vh->info().default_z)
          return false;
        if (vh == end) break;

        face.push_back(Point_3(q.x(), q.y(), vh->info().z));
        if (num_iters == m_max_num_iters)
          return false;
        ++num_iters;
      } while (vh != end);
      return is_valid_traversal(face);
    }

    bool find_first_edge(
      const Face_handle& fh,
      const CDT& cdt,
      Edge& edge) const {

      for (std::size_t k = 0; k < 3; ++k) {
        edge = std::make_pair(fh, k);
        if (cdt.is_constrained(edge))
          return true;
      }
      return false;
    }

    void get_next_vertex_handle(
      const CDT& cdt,
      Vertex_handle& vh,
      Edge& edge) const {

      const std::size_t idx = edge.first->index(vh);
      Edge next = std::make_pair(edge.first, (idx + 2) % 3);
      while (!cdt.is_constrained(next)) {

        const auto fh = next.first->neighbor(next.second);
        const auto tmp = next.first->vertex((next.second + 1) % 3);
        const std::size_t tmp_idx = fh->index(tmp);
        next = std::make_pair(fh, (tmp_idx + 2) % 3);
      }
      vh = next.first->vertex((next.second + 2) % 3);
      edge = next;
    }

    bool is_valid_traversal(
      const Polygon& face) const {

      if (face.size() < 3) return false;
      for (std::size_t i = 0; i < face.size(); ++i) {
        const auto& p = face[i];
        for (std::size_t j = 0; j < face.size(); ++j) {
          const auto& q = face[j];

          if (i == j) continue;
          if (internal::are_equal_points_3(p, q))
            return false;
        }
      }
      return true;
    }

    void fix_orientation(
      Polygon& face) const {

      std::vector<Point_2> polygon_2;
      polygon_2.reserve(face.size());
      for (const auto& p : face)
        polygon_2.push_back(Point_2(p.x(), p.y()));
      if (CGAL::orientation_2(
      polygon_2.begin(), polygon_2.end()) == CGAL::CLOCKWISE)
        std::reverse(face.begin(), face.end());
    }
  };

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_INTERNAL_WALLS_ESTIMATOR_H
