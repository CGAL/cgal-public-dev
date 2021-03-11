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

#ifndef CGAL_LEVELS_OF_DETAIL_INTERNAL_WALLS_ESTIMATOR_EXP_1_H
#define CGAL_LEVELS_OF_DETAIL_INTERNAL_WALLS_ESTIMATOR_EXP_1_H

// Example:

// // Add extra walls, which originate from the roof boundaries.
// westimator.add_extra_walls(
//   m_cluster,
//   m_data.point_map_3,
//   m_roof_points_3,
//   m_data.parameters.buildings.alpha_shape_size_2,
//   m_building_walls);

// STL includes.
#include <vector>
#include <utility>

// CGAL includes.
#include <CGAL/assertions.h>
#include <CGAL/Optimal_transportation_reconstruction_2.h>

// Internal includes.
#include <CGAL/Levels_of_detail/internal/utils.h>
#include <CGAL/Levels_of_detail/internal/struct.h>

// Spatial search.
#include <CGAL/Levels_of_detail/internal/Spatial_search/Nearest_face_neighbor_query.h>

// Shape detection.
#include <CGAL/Levels_of_detail/internal/Shape_detection/Region_growing.h>
#include <CGAL/Levels_of_detail/internal/Shape_detection/Coplanar_region.h>

// Simplification.
#include <CGAL/Levels_of_detail/internal/Simplification/Alpha_shapes_filtering_2.h>

// Testing.
#include "../../../../../test/Levels_of_detail/include/Saver.h"

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

  template<typename GeomTraits>
  class Building_walls_estimator_exp_1 {

  public:
    using Traits = GeomTraits;

    using FT = typename Traits::FT;
    using Point_2 = typename Traits::Point_2;
    using Point_3 = typename Traits::Point_3;
    using Vector_3 = typename Traits::Vector_3;
    using Plane_3 = typename Traits::Plane_3;
    using Segment_2 = typename Traits::Segment_2;

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

    Building_walls_estimator_exp_1(
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

    // These extra walls originate from the roof boundaries.
    template<
    typename InputRange,
    typename PointMap3>
    void add_extra_walls(
      const InputRange& input_range,
      const PointMap3 point_map_3,
      const std::vector<Indices>& roof_points_3,
      const FT alpha_shape_size_2,
      std::vector<Approximate_face>& walls) const {

      using Color = CGAL::Color;
      using Saver = Saver<Traits>;
      using Identity_map_3 = CGAL::Identity_property_map<Point_3>;

      std::vector<Point_2> points2;
      std::vector<Point_3> points3;
      std::vector<Segment_2> segments;

      std::size_t psi = 0;
      for (const auto& region : roof_points_3) {
        points3.clear();
        for (const std::size_t idx : region) {
          const auto& p = get(point_map_3, *(input_range.begin() + idx));
          points3.push_back(p);
        }

        // Point_3 b;
        // internal::compute_barycenter_3(points3, b);
        // Plane_3 plane; Identity_map_3 identity_map_3;
        // internal::plane_from_points_3(points3, identity_map_3, plane);

        points2.clear();
        for (const auto& p : points3)
          points2.push_back(internal::point_2_from_point_3(p));
        apply_filtering(alpha_shape_size_2, points2);

        segments.clear();
        const FT scale = alpha_shape_size_2 / FT(2);
        compute_optimal_transport(
          points2,
          scale,
          segments);

        Approximate_face face;
        for (const auto& segment : segments) {
          const auto& s = segment.source();
          const auto& t = segment.target();

          const Point_3 p1 = Point_3(s.x(), s.y(), m_bottom_z);
          const Point_3 p2 = Point_3(t.x(), t.y(), m_bottom_z);
          const Point_3 p3 = Point_3(t.x(), t.y(), m_top_z);
          const Point_3 p4 = Point_3(s.x(), s.y(), m_top_z);

          face.polygon.clear();
          face.polygon.reserve(4);

          face.polygon.push_back(p1);
          face.polygon.push_back(p2);
          face.polygon.push_back(p3);
          face.polygon.push_back(p4);

          walls.push_back(face);
        }

        // Saver saver;
        // std::vector<Point_3> tmp;
        // for (const auto& p : points2)
        //   tmp.push_back(Point_3(p.x(), p.y(), FT(0)));
        // const Color color(0, 0, 0);
        // const std::string pname =
        // "/Users/monet/Documents/lod/logs/buildings/groups/point-set-" + std::to_string(psi);
        // saver.export_points(tmp, color, pname);

        // std::vector< std::vector<Point_3> > polylines(segments.size());
        // for (std::size_t i = 0; i < segments.size(); ++i) {
        //   const auto& s = segments[i].source();
        //   const auto& t = segments[i].target();
        //   polylines[i].push_back(Point_3(s.x(), s.y(), FT(0)));
        //   polylines[i].push_back(Point_3(t.x(), t.y(), FT(0)));
        // }

        // const std::string lname =
        // "/Users/monet/Documents/lod/logs/buildings/polylines/polylines-" + std::to_string(psi);
        // saver.export_polylines(polylines, lname);
        // ++psi;
      }
    }

  private:
    const std::vector<Boundary>& m_boundaries;
    const FT m_bottom_z;
    const FT m_top_z;
    const std::size_t m_max_num_iters;

    void apply_filtering(
      const FT alpha_shape_size_2,
      std::vector<Point_2>& points2) const {

      using Identity_map_2 = CGAL::Identity_property_map<Point_2>;
      using Alpha_shapes_filtering_2 = internal::Alpha_shapes_filtering_2<Traits>;

      const std::size_t nump = points2.size();
      Alpha_shapes_filtering_2 filtering(alpha_shape_size_2);
      const FT sampling_2 = alpha_shape_size_2 / FT(2);

      Identity_map_2 identity_map_2;
      filtering.add_points(points2, identity_map_2);
      points2.clear();
      filtering.get_filtered_points(sampling_2, points2);
    }

    void compute_optimal_transport(
      const std::vector<Point_2>& points2,
      const FT scale,
      std::vector<Segment_2>& segments) const {

      using Identity_map_2 = CGAL::Identity_property_map<Point_2>;
      using Otr = CGAL::Optimal_transportation_reconstruction_2<Traits, Identity_map_2>;

      Identity_map_2 identity_map_2;
      Otr otr(points2, identity_map_2);
      otr.run_under_wasserstein_tolerance(scale);
      otr.list_output(
        boost::make_function_output_iterator([&](const Point_2&) -> void { }),
        boost::make_function_output_iterator([&](const Segment_2& segment_2) -> void {
          segments.push_back(segment_2);
        })
      );
    }

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

      for (const auto& region : regions)
        merge_faces(faces, region, walls);
    }

    void merge_faces(
      const std::vector<Polygon>& faces,
      const Indices& region,
      std::vector<Approximate_face>& merged) const {

      CGAL_assertion(region.size() > 0);

      Vector_3 m;
      bool success = internal::compute_normal_3(faces[region[0]], m);
      if (!success) return;
      const Vector_3 n = Vector_3(FT(0), FT(0), FT(1));
      if (m == -n) m = n;

      FT angle_3d; Vector_3 axis;
      success = internal::compute_angle_and_axis_3(
        m, n, angle_3d, axis);
      if (!success) return;
      const FT angle_deg = angle_3d * FT(180) / static_cast<FT>(CGAL_PI);

      Point_3 b;
      internal::compute_barycenter_3(faces, region, b);

      std::vector<Polygon> rotated;
      if (angle_deg != FT(0) && angle_deg != FT(180))
        internal::rotate_polygons_3(faces, region, angle_3d, axis, b, rotated);

      CDT cdt;
      triangulate(rotated, cdt);
      if (cdt.number_of_faces() == 0)
        return;

      Polygon merged_face;
      success = create_merged_face(cdt, merged_face);
      if (!success || merged_face.size() < 3) return;
      fix_orientation(merged_face);

      Approximate_face face;
      if (angle_deg != FT(0) && angle_deg != FT(180))
        internal::rotate_polygon_3(merged_face, -angle_3d, axis, b, face.polygon);
      merged.push_back(face);
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

#endif // CGAL_LEVELS_OF_DETAIL_INTERNAL_WALLS_ESTIMATOR_EXP_1_H
