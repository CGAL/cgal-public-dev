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

#ifndef CGAL_LEVELS_OF_DETAIL_INTERNAL_SMOOTH_GROUND_BUILDER_H
#define CGAL_LEVELS_OF_DETAIL_INTERNAL_SMOOTH_GROUND_BUILDER_H

#include <CGAL/license/Levels_of_detail.h>

// STL includes.
#include <set>
#include <map>
#include <vector>
#include <memory>

// CGAL includes.
#include <CGAL/barycenter.h>
#include <CGAL/Kernel/global_functions.h>
#include <CGAL/assertions.h>

// Internal includes.
#include <CGAL/Levels_of_detail/internal/utils.h>
#include <CGAL/Levels_of_detail/internal/Ground/Planar_ground_builder.h>

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

template<
  typename GroundBase,
  typename PointRange,
  typename PointMap,
  typename NeighborQuery>
  class Smooth_ground_builder : public Planar_ground_builder<GroundBase> {

  public:
    using Ground_base = GroundBase;
    using Point_range = PointRange;
    using Point_map = PointMap;
    using Neighbor_query = NeighborQuery;

    using Traits = typename Ground_base::Traits;
    using FT = typename Traits::FT;
    using Point_2 = typename Traits::Point_2;
    using Point_3 = typename Traits::Point_3;
    using Plane_3 = typename Traits::Plane_3;
    using Triangle_3 = typename Traits::Triangle_3;

    using Triangulation = typename Ground_base::Triangulation::Delaunay;
    using Face_handle = typename Triangulation::Face_handle;
    using Vertex_handle = typename Triangulation::Vertex_handle;
    using Point_iterator = typename Point_range::const_iterator;
    using m_base = Planar_ground_builder<Ground_base>;
    using Location_type = typename Triangulation::Locate_type;

    struct Candidate_face {
      Face_handle face;
      FT max_error;
      std::vector<Point_iterator> inliers;
      Candidate_face(Face_handle fh = Face_handle()) :
      face(fh), max_error(FT(0)) { }
    };

    using Candidate_face_ptr = std::shared_ptr<Candidate_face>;
    struct Compare_candidates {
      bool operator()(
        const Candidate_face_ptr& a, const Candidate_face_ptr& b) const {
        if (a->max_error != b->max_error)
          return a->max_error > b->max_error;
        return a->face < b->face;
      }
    };

    using Face_map = std::map<Face_handle, Candidate_face_ptr>;
    using Face_queue = std::set<Candidate_face_ptr, Compare_candidates>;

    Smooth_ground_builder(
      Ground_base& ground_base,
      const Point_range& points,
      const Point_map point_map,
      const Neighbor_query& neighbor_query,
      const FT ground_precision,
      const bool scale_bbox = false) :
    m_base(ground_base),
    m_points(points),
    m_point_map(point_map),
    m_neighbor_query(neighbor_query),
    m_ground_precision(ground_precision),
    m_scale_bbox(scale_bbox)
    { }

    void finilize() {

      set_real_heights();
      m_base::set_face_heights();
      refine();
      if (m_scale_bbox)
        scale_bbox(FT(5) / FT(4));
    }

  private:
    const Point_range& m_points;
    const Point_map m_point_map;
    const Neighbor_query& m_neighbor_query;
    const FT m_ground_precision;
    const bool m_scale_bbox;

    void scale_bbox(const FT scale) {

      Triangulation& tri = m_base::m_ground_base.triangulation.delaunay;
      auto bbox = m_base::m_ground_base.bbox;
      internal::scale_polygon_2(scale, bbox);

      for (std::size_t i = 0; i < bbox.size(); ++i) {
        const std::size_t ip = (i + 1) % bbox.size();

        const Point_2& p = bbox[i];
        const Point_2& q = bbox[ip];

        const Point_2 c = internal::middle_point_2(p, q);
        tri.insert(p); tri.insert(c);
      }

      std::vector<std::size_t> neighbors;
      for (auto vh = tri.finite_vertices_begin();
      vh != tri.finite_vertices_end(); ++vh)
        if (vh->info().z == vh->info().default_z)
          vh->info().z = get_z(vh, neighbors);
      m_base::set_face_heights();
    }

    void set_real_heights() {

      Triangulation& tri = m_base::m_ground_base.triangulation.delaunay;
      std::vector<std::size_t> neighbors;
      for (auto vh = tri.finite_vertices_begin();
      vh != tri.finite_vertices_end(); ++vh)
        vh->info().z = get_z(vh, neighbors);
    }

    template<typename FH>
    void set_ground_heights(FH& fh) const {

      std::vector<std::size_t> neighbors;
      for (std::size_t k = 0; k < 3; ++k) {

        const Vertex_handle& vh = fh->vertex(k);
        const FT z = get_z(vh, neighbors);
        const std::size_t idx = fh->index(vh);
        CGAL_assertion(idx >= 0 && idx < 3);
        fh->info().z[idx] = z;
      }
    }

    FT get_z(
      const Vertex_handle& vh,
      std::vector<std::size_t>& neighbors) const {

      m_neighbor_query(vh->point(), neighbors);
      FT z = FT(0);
      for (const std::size_t nidx : neighbors)
        z += get(m_point_map, *(m_points.begin() + nidx)).z();
      z /= static_cast<FT>(neighbors.size());
      return z;
    }

    void refine() {
      const FT tolerance = m_ground_precision;
      Face_map face_map; Face_queue todo;
      initialize_queue(face_map, todo, m_ground_precision);
      auto& tri = m_base::m_ground_base.triangulation.delaunay;

      Face_handle hint;
      while (!todo.empty()) {
        const Candidate_face_ptr candidate = *(todo.begin());
        todo.erase(todo.begin());

        const bool out_of_tolerance = (candidate->max_error > tolerance * tolerance);
        const bool badly_shaped = !well_shaped(candidate->face);

        if (!out_of_tolerance && !badly_shaped) continue;
        if (is_too_small(candidate->face, FT(3) * tolerance)) continue;

        // Get circumcenter and conflict zone.
        Point_2 center = CGAL::circumcenter(
          candidate->face->vertex(0)->point(),
          candidate->face->vertex(1)->point(),
          candidate->face->vertex(2)->point());

        Location_type type; int stub;
        hint = tri.locate(center, type, stub, hint);
        if (type == Triangulation::VERTEX || !is_valid_face(tri, hint)) {
          if (out_of_tolerance) {
            center = CGAL::barycenter(
              candidate->face->vertex(0)->point(), FT(1),
              candidate->face->vertex(1)->point(), FT(1),
              candidate->face->vertex(2)->point(), FT(1));

            hint = tri.locate(center, type, stub, hint);
            if (type == Triangulation::VERTEX || !is_valid_face(tri, hint))
            continue;
          } else continue;
        }
        std::vector<Face_handle> conflict;
        tri.get_conflicts(center, std::back_inserter(conflict));

        // Recover points and remove conflict cells from local structures.
        std::vector<Point_iterator> points;
        for (std::size_t i = 0; i < conflict.size(); ++i) {
          if (!is_valid_face(tri, conflict[i])) continue;

          const auto filter = face_map.find(conflict[i]);
          if (filter == face_map.end()) continue;

          const Candidate_face_ptr cface = filter->second;
          if (!(cface->inliers.empty()))
            std::copy(cface->inliers.begin(), cface->inliers.end(),
              std::back_inserter(points));

          face_map.erase(filter);
          todo.erase(cface);
        }

        // Insert new vertex.
        const Vertex_handle v = tri.insert(center, hint);
        std::vector<Candidate_face_ptr> new_faces;

        auto circ = tri.incident_faces(v);
        const auto start = circ;
        do {
          circ->info().urban_tag = Urban_object_type::GROUND;
          set_ground_heights(circ);
          const Candidate_face_ptr cface = std::make_shared<Candidate_face>(circ);
          face_map.insert(std::make_pair(circ, cface));
          new_faces.push_back(cface);
          ++circ;
        }
        while (circ != start);

        // Redistribute points.
        for (auto pit = points.begin(); pit != points.end(); ++pit) {
          const Point_3& point_3 = get(m_point_map, **pit);
          const Point_2 point_2 = internal::point_2_from_point_3(point_3);

          hint = tri.locate(point_2, hint);
          const auto filter = face_map.find(hint);
          CGAL_assertion(filter != face_map.end());
          Candidate_face_ptr cface = filter->second;

          Triangle_3 triangle;
          internal::triangle_3(hint, triangle);
          const FT sq_dist = CGAL::squared_distance(point_3, triangle);

          cface->inliers.push_back(*pit);
          cface->max_error = CGAL::max(cface->max_error, sq_dist);
        }

        // Insert new faces.
        for (std::size_t i = 0; i < new_faces.size(); ++i)
          todo.insert(new_faces[i]);
      }
    }

    void initialize_queue(
      Face_map& face_map, Face_queue& todo, FT tolerance) {
      auto& tri = m_base::m_ground_base.triangulation.delaunay;

      for(auto fit = tri.finite_faces_begin();
      fit != tri.finite_faces_end(); ++fit) {
        if (!is_valid_face(tri, fit)) continue;
        face_map.insert(std::make_pair(
          fit, std::make_shared<Candidate_face>(fit)));
      }

      Face_handle hint;
      for (auto pit = m_points.begin(); pit != m_points.end(); ++pit) {
        const Point_3& point_3 = get(m_point_map, *pit);
        const Point_2 point_2 = internal::point_2_from_point_3(point_3);

        hint = tri.locate(point_2, hint);
        if (!is_valid_face(tri, hint)) continue;

        const auto filter = face_map.find(hint);
        CGAL_assertion(filter != face_map.end());
        Candidate_face_ptr candidate = filter->second;

        Triangle_3 triangle;
        internal::triangle_3(hint, triangle);
        const FT sq_dist = CGAL::squared_distance(point_3, triangle);

        candidate->inliers.push_back(pit);
        candidate->max_error = CGAL::max(candidate->max_error, sq_dist);
      }

      for (auto fit = face_map.begin(); fit != face_map.end(); ++fit)
        todo.insert(fit->second);
    }

    bool well_shaped(const Face_handle& fh) const {
      const Point_2& pa = fh->vertex(0)->point();
      const Point_2& pb = fh->vertex(1)->point();
      const Point_2& pc = fh->vertex(2)->point();

      double area =
      2.0 * CGAL::to_double(CGAL::area(pa, pb, pc));area = area * area;
      const double a = CGAL::to_double(CGAL::squared_distance(pb, pc));
      const double b = CGAL::to_double(CGAL::squared_distance(pc, pa));
      const double c = CGAL::to_double(CGAL::squared_distance(pa, pb));

      double q;
      if(a < b) {
        if(a < c)
          q = area / (b * c);
        else
          q = area / (a * b);
      } else {
        if(b < c)
          q = area / (a * c);
        else
          q = area / (a * b);
      }
      return q > 0.125;
    }

    bool is_too_small(
      const Face_handle& fh,
      const FT min_size) const {
      for (std::size_t ka = 0; ka < 3; ++ka) {
        const std::size_t kb = (ka + 1) % 3;
        const Point_2& a = fh->vertex(ka)->point();
        const Point_2& b = fh->vertex(kb)->point();
        if (CGAL::squared_distance(a, b) > min_size * min_size)
          return false;
      }
      return true;
    }

    bool is_valid_face(
      const Triangulation& tri, const Face_handle& fh) const {
      return !tri.is_infinite(fh) &&
      fh->info().urban_tag == Urban_object_type::GROUND;
    }
  };

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_INTERNAL_SMOOTH_GROUND_BUILDER_H
