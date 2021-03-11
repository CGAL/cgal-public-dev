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

#ifndef CGAL_LEVELS_OF_DETAIL_INTERNAL_WALLS_ESTIMATOR_EXP_2_H
#define CGAL_LEVELS_OF_DETAIL_INTERNAL_WALLS_ESTIMATOR_EXP_2_H

// Example:

// westimator.add_extra_walls(
//   m_cluster,
//   m_data.point_map_3,
//   m_roof_points_3,
//   m_data.parameters.buildings.region_growing_min_area_3,
//   m_data.parameters.noise_level,
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
#include <CGAL/Levels_of_detail/internal/Spatial_search/K_neighbor_query.h>
#include <CGAL/Levels_of_detail/internal/Spatial_search/Nearest_face_neighbor_query.h>

// Shape detection.
#include <CGAL/Levels_of_detail/internal/Shape_detection/Region_growing.h>
#include <CGAL/Levels_of_detail/internal/Shape_detection/Coplanar_region.h>

// Testing.
#include "../../../../../test/Levels_of_detail/include/Saver.h"

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

  template<typename GeomTraits>
  class Building_walls_estimator_exp_2 {

  public:
    using Traits = GeomTraits;

    using FT = typename Traits::FT;
    using Point_2 = typename Traits::Point_2;
    using Point_3 = typename Traits::Point_3;
    using Vector_3 = typename Traits::Vector_3;
    using Segment_2 = typename Traits::Segment_2;
    using Triangle_2 = typename Traits::Triangle_2;

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
    using Face_iterator = typename CDT::Finite_faces_iterator;
    using Face_handle = typename CDT::Face_handle;
    using Vertex_iterator = typename CDT::Finite_vertices_iterator;
    using Vertex_handle = typename CDT::Vertex_handle;
    using Vertex_circulator = typename CDT::Vertex_circulator;
    using Edge = typename CDT::Edge;

    using Vhs = std::vector<Vertex_handle>;
    using Vh_pair = std::pair<Vertex_handle, Vertex_handle>;

    using Pair_2 = std::pair<Point_2, std::pair<std::size_t, FT> >;
    using Pairs_2 = std::vector<Pair_2>;

    using Point_map_2 = CGAL::First_of_pair_property_map<Pair_2>;
    using KNQ = internal::K_neighbor_query<Traits, Pairs_2, Point_map_2>;

    using Color = CGAL::Color;
    using Saver = Saver<Traits>;

    Building_walls_estimator_exp_2(
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
      const FT min_height,
      const FT noise,
      std::vector<Approximate_face>& walls) const {

      std::cout << "height: " << min_height << " and noise: " << noise << std::endl;

      // Create points.
      Pairs_2 points2;
      create_points2(input_range, point_map_3, roof_points_3, points2);
      std::vector<Point_2> pts2;

      // Version 1.
      // CDT cdt;
      // create_triangulation(points2, cdt);
      // extract_boundary_points_from_cdt_vertices(cdt, min_height, pts2);

      // Version 2.
      // CDT cdt;
      // create_triangulation(points2, cdt);
      // extract_boundary_points_from_cdt_faces(cdt, min_height, pts2);

      // Version 3.
      Point_map_2 point_map_2;
      KNQ knq(points2, FT(3), point_map_2);
      extract_boundary_points_using_kd_tree(points2, knq, min_height, pts2);

      // Create segments.
      std::vector<Segment_2> segments;
      compute_optimal_transport(pts2, noise, segments);

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

      // Saving.
     Saver saver;

      // Triangulation.
      // std::vector< std::vector<Point_3> > trisave;
      // for (auto fit = cdt.finite_faces_begin();
      // fit != cdt.finite_faces_end(); ++fit) {

      //   const auto& a = fit->vertex(0)->point();
      //   const auto& b = fit->vertex(1)->point();
      //   const auto& c = fit->vertex(2)->point();

      //   const Point_3 p1 = Point_3(a.x(), a.y(), FT(0));
      //   const Point_3 p2 = Point_3(b.x(), b.y(), FT(0));
      //   const Point_3 p3 = Point_3(c.x(), c.y(), FT(0));

      //   std::vector<Point_3> poly;
      //   poly.push_back(p1);
      //   poly.push_back(p2);
      //   poly.push_back(p3);
      //   trisave.push_back(poly);
      // }
      // const Color color1(255, 0, 0);
      // const std::string tname =
      // "/Users/monet/Documents/lod/logs/buildings/triangulation";
      // saver.export_polygon_soup(trisave, color1, tname);

      // Points.
      std::vector<Point_3> tmp;
      for (const auto& p : pts2)
        tmp.push_back(Point_3(p.x(), p.y(), FT(0)));
      const Color color2(0, 0, 0);
      const std::string pname =
      "/Users/monet/Documents/lod/logs/buildings/points";
      saver.export_points(tmp, color2, pname);

      // Polylines.
      std::vector< std::vector<Point_3> > polylines(segments.size());
      for (std::size_t i = 0; i < segments.size(); ++i) {
        const auto& s = segments[i].source();
        const auto& t = segments[i].target();
        polylines[i].push_back(Point_3(s.x(), s.y(), FT(0)));
        polylines[i].push_back(Point_3(t.x(), t.y(), FT(0)));
      }
      const std::string lname =
      "/Users/monet/Documents/lod/logs/buildings/segments";
      saver.export_polylines(polylines, lname);
    }

  private:
    const std::vector<Boundary>& m_boundaries;
    const FT m_bottom_z;
    const FT m_top_z;
    const std::size_t m_max_num_iters;

    template<
    typename InputRange,
    typename PointMap3>
    void create_points2(
      const InputRange& input_range,
      const PointMap3 point_map_3,
      const std::vector<Indices>& roof_points_3,
      Pairs_2& points2) const {

      points2.clear();
      std::size_t nump = 0;
      for (const auto& region : roof_points_3)
        nump += region.size();
      points2.reserve(nump);

      for (std::size_t i = 0; i < roof_points_3.size(); ++i) {
        for (const std::size_t idx : roof_points_3[i]) {
          const auto& p = get(point_map_3, *(input_range.begin() + idx));
          points2.push_back(
            std::make_pair(
              Point_2(p.x(), p.y()),
              std::make_pair(i, p.z())));
        }
      }
    }

    void create_triangulation(
      const Pairs_2& points2,
      CDT& cdt) const {

      cdt.clear();
      for (const auto& pair : points2) {
        auto vh = cdt.insert(pair.first);
        vh->info().object_index = pair.second.first;
        vh->info().z = pair.second.second;
      }
    }

    void extract_boundary_points_from_cdt_vertices(
      const CDT& cdt,
      const FT min_height,
      std::vector<Point_2>& pts2) const {

      std::vector<Vertex_circulator> vhs;
      for (auto vh1 = cdt.finite_vertices_begin();
      vh1 != cdt.finite_vertices_end(); ++vh1) {
        Vertex_circulator vc = cdt.incident_vertices(vh1);
        get_vhs(vc, vhs);
        extract_boundary_point_from_cdt_vertices(cdt, vh1, vhs, min_height, pts2);
      }
    }

    void extract_boundary_point_from_cdt_vertices(
      const CDT& cdt,
      const Vertex_iterator& curr,
      const std::vector<Vertex_circulator>& vhs,
      const FT min_height,
      std::vector<Point_2>& pts2) const {

      for (std::size_t i = 0; i < vhs.size(); ++i) {
        if (should_be_added(cdt, curr, vhs[i], min_height)) {
          pts2.push_back(curr->point());
          return;
        }
      }
    }

    void get_vhs(
      Vertex_circulator& vc,
      std::vector<Vertex_circulator>& vhs) const {

      vhs.clear();
      if (vc.is_empty()) return;
      Vertex_circulator end = vc;
      do {
        vhs.push_back(vc); ++vc;
      } while (vc != end);
    }

    bool should_be_added(
      const CDT& cdt,
      const Vertex_iterator& vh1,
      const Vertex_circulator& vh2,
      const FT min_height) const {

      const bool cond1 = apply_conditions(
        vh1->info().object_index, vh2->info().object_index,
        vh1->info().z, vh2->info().z,
        min_height);
      const bool cond2 = ( !cdt.is_infinite(vh1) && !cdt.is_infinite(vh2) );

      return (cond1 && cond2);
    }

    bool apply_conditions(
      const std::size_t index1, const std::size_t index2,
      const FT z1, const FT z2,
      const FT min_height) const {

      const bool cond1 = ( index1 != std::size_t(-1) && index2 != std::size_t(-1) );
      const bool cond2 = ( index1 != index2 );
      const bool cond3 = ( CGAL::abs(z1 - z2) > min_height );

      return (cond1 && cond2 && cond3);
    }

    void extract_boundary_points_using_kd_tree(
      const Pairs_2& points2, KNQ& knq, const FT min_height,
      std::vector<Point_2>& pts2) const {

      Indices neighbors;
      for (const auto& p : points2) {
        knq(p.first, neighbors);
        for (const std::size_t idx : neighbors) {
          const auto& q = points2[idx];
          if (apply_conditions(
            p.second.first, q.second.first,
            p.second.second, q.second.second,
            min_height)) {

            pts2.push_back(p.first);
            break;
          }
        }
      }
    }

    void extract_boundary_points_from_cdt_faces(
      const CDT& cdt,
      const FT min_height,
      std::vector<Point_2>& pts2) const {

      tag_faces(cdt, min_height);
      clean_faces(cdt);

      create_points(cdt, pts2);
    }

    void tag_faces(
      const CDT& cdt,
      const FT min_height) const {

      for (auto fh = cdt.finite_faces_begin();
      fh != cdt.finite_faces_end(); ++fh) {
        if (is_necessary_face(fh, min_height))
          fh->info().tagged = true;
        else
          fh->info().tagged = false;
      }
    }

    bool is_necessary_face(
      const Face_iterator& fh,
      const FT min_height) const {

      const auto& vh1 = fh->vertex(0);
      const auto& vh2 = fh->vertex(1);
      const auto& vh3 = fh->vertex(2);
      return
        face_should_be_added(vh1, vh2, min_height) ||
        face_should_be_added(vh2, vh3, min_height) ||
        face_should_be_added(vh3, vh1, min_height);
    }

    bool face_should_be_added(
      const Vertex_handle& vh1,
      const Vertex_handle& vh2,
      const FT min_height) const {

      return apply_conditions(
        vh1->info().object_index, vh2->info().object_index,
        vh1->info().z, vh2->info().z,
        min_height);
    }

    void clean_faces(
      const CDT& cdt) const {

      close_gaps(cdt);
      close_single_faces(cdt);
    }

    void close_gaps(
      const CDT& cdt) const {

      for (auto fh = cdt.finite_faces_begin();
      fh != cdt.finite_faces_end(); ++fh) {
        if (is_gap_face(fh))
          fh->info().tagged = true;
      }
    }

    bool is_gap_face(
      const Face_iterator& fh) const {

      const auto& fh1 = fh->neighbor(0);
      const auto& fh2 = fh->neighbor(1);
      const auto& fh3 = fh->neighbor(2);

      return  (fh1->info().tagged && fh2->info().tagged) ||
              (fh2->info().tagged && fh3->info().tagged) ||
              (fh3->info().tagged && fh1->info().tagged);
    }

    void close_single_faces(
      const CDT& cdt) const {

      for (auto fh = cdt.finite_faces_begin();
      fh != cdt.finite_faces_end(); ++fh) {
        if (is_single_face(fh))
          fh->info().tagged = false;
      }
    }

    bool is_single_face(
      const Face_iterator& fh) const {

      const auto& fh1 = fh->neighbor(0);
      const auto& fh2 = fh->neighbor(1);
      const auto& fh3 = fh->neighbor(2);

      return  (!fh1->info().tagged && !fh2->info().tagged) ||
              (!fh2->info().tagged && !fh3->info().tagged) ||
              (!fh3->info().tagged && !fh1->info().tagged);
    }

    void create_points(
      const CDT& cdt,
      std::vector<Point_2>& pts2) const {

      for (auto fh = cdt.finite_faces_begin();
      fh != cdt.finite_faces_end(); ++fh) {
        if (is_valid_face(fh)) {

          const auto& vh1 = fh->vertex(0);
          const auto& vh2 = fh->vertex(1);
          const auto& vh3 = fh->vertex(2);

          pts2.push_back(vh1->point());
          pts2.push_back(vh2->point());
          pts2.push_back(vh3->point());
        }
      }
    }

    bool is_valid_face(
      const Face_iterator& fh) const {

      const Point_2& p1 = fh->vertex(0)->point();
      const Point_2& p2 = fh->vertex(1)->point();
      const Point_2& p3 = fh->vertex(2)->point();

      const Triangle_2 triangle = Triangle_2(p1, p2, p3);
      return ( fh->info().tagged && triangle.area() < (FT(1) / FT(10)) );
    }

    void compute_optimal_transport(
      const std::vector<Point_2>& pts2,
      const FT noise,
      std::vector<Segment_2>& segments) const {

      using PMap_2 = CGAL::Identity_property_map<Point_2>;
      using Otr = CGAL::Optimal_transportation_reconstruction_2<Traits, PMap_2>;

      PMap_2 pmap2;
      Otr otr(pts2, pmap2);
      otr.run_under_wasserstein_tolerance(noise * FT(2));
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

#endif // CGAL_LEVELS_OF_DETAIL_INTERNAL_WALLS_ESTIMATOR_EXP_2_H
