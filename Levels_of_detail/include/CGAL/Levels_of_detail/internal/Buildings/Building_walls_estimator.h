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
// #include <CGAL/number_utils.h>
// #include <CGAL/Polygon_2_algorithms.h>
// #include <CGAL/Triangulation_face_base_with_info_2.h>
// #include <CGAL/Constrained_Delaunay_triangulation_2.h>
// #include <CGAL/Constrained_triangulation_face_base_2.h>
// #include <CGAL/Triangulation_vertex_base_with_info_2.h>

// Internal includes.
#include <CGAL/Levels_of_detail/internal/utils.h>
#include <CGAL/Levels_of_detail/internal/struct.h>

// Spacial search.
#include <CGAL/Levels_of_detail/internal/Spacial_search/Nearest_face_neighbor_query.h>

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

    /*
    using Segment_2 = typename Traits::Segment_2;
    */

    using Indices = std::vector<std::size_t>;
    using Boundary = internal::Boundary<Traits>;
    using Approximate_face = internal::Partition_edge_3<Traits>;
    using Polygon = std::vector<Point_3>;

    using Nearest_face_neighbor_query = internal::Nearest_face_neighbor_query<Traits>;
    using Coplanar_region = internal::Coplanar_region<Traits>;
    using Region_growing = internal::Region_growing<
    std::vector<Polygon>, Nearest_face_neighbor_query, Coplanar_region>;

    /*
    using Face_info = Face_info<Traits>;
    using Vertex_info = Vertex_info<Traits>;

    using VB = CGAL::Triangulation_vertex_base_with_info_2<Vertex_info, Traits>;
    using FB = CGAL::Triangulation_face_base_with_info_2<Face_info, Traits>;
    
    using CFB = CGAL::Constrained_triangulation_face_base_2<Traits, FB>;
    using TAG = CGAL::Exact_predicates_tag;
    using TDS = CGAL::Triangulation_data_structure_2<VB, CFB>;

    using CDT = CGAL::Constrained_Delaunay_triangulation_2<Traits, TDS, TAG>;
    
    using Face_handle = typename CDT::Face_handle;
    using Vertex_handle = typename CDT::Vertex_handle;
    using Edge = typename CDT::Edge; 
    */

    Building_walls_estimator(
      const std::vector<Boundary>& boundaries,
      const FT bottom_z,
      const FT top_z) :
    m_boundaries(boundaries),
    m_bottom_z(bottom_z),
    m_top_z(top_z)
    /*,
    m_tolerance(FT(1) / FT(100000)),
    m_default_height(-FT(100000000000000)),
    m_max_num_iters(100) */
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

  private:
    const std::vector<Boundary>& m_boundaries;
    const FT m_bottom_z;
    const FT m_top_z;

    /*
    const FT m_tolerance;
    const FT m_default_height;
    const std::size_t m_max_num_iters; */

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
      walls.reserve(regions.size());

      for (const auto& region : regions)
        merge_faces(faces, region, walls);
    }

    void merge_faces(
      const std::vector<Polygon>& faces,
      const Indices& region,
      std::vector<Approximate_face>& merged) const {

      Vector_3 source_normal;
      bool success = compute_normal(faces, region, source_normal);
      if (!success) {
        std::cerr << "Error: source normal!" << std::endl;
        exit(EXIT_FAILURE); }
      const Vector_3 target_normal = Vector_3(FT(0), FT(0), FT(1));

      if (source_normal == -target_normal) 
        source_normal = target_normal;

      FT angle; Vector_3 axis;
      success = compute_angle_and_axis(
        source_normal, target_normal, angle, axis);
      if (!success) {
        std::cerr << "Error: angle and axis!" << std::endl;   
        exit(EXIT_FAILURE); }

      Point_3 b;
      compute_barycenter(faces, region, b);
      const FT angle_deg = angle * FT(180) / static_cast<FT>(CGAL_PI);
                
      /*
      std::vector< std::vector<Point_3> > rotated;
      if (angle_deg != FT(0) && angle_deg != FT(180))
        rotate_walls(walls, region, angle, axis, b, rotated);
      
      CDT cdt;
      triangulate_region_facets(rotated, cdt);

      std::vector<Point_3> final_face;
      if (cdt.number_of_faces() != 0) 
        get_back_region_facets(cdt, final_face);

      fix_orientation(final_face);

      std::vector< std::vector<Point_3> > region_facets = { final_face };
      std::vector<std::size_t> indices = { 0 };

      rotated.clear();
      if (angle_deg != FT(0) && angle_deg != FT(180))
        rotate_walls(region_facets, indices, -angle, axis, b, rotated);

      result.push_back(rotated[0]); */
    }

    /*
    void rotate_walls(
      const std::vector< std::vector<Point_3> >& walls,
      const std::vector<std::size_t>& region, 
      const FT angle, 
      const Vector_3& axis, 
      const Point_3& b,
      std::vector< std::vector<Point_3> >& rotated) const {

      rotated.resize(region.size());
      for (std::size_t i = 0; i < region.size(); ++i)
        rotate_wall(walls[region[i]], angle, axis, b, rotated[i]);
    }

    void rotate_wall(
      const std::vector<Point_3>& wall, 
      const FT angle, 
      const Vector_3& axis, 
      const Point_3& b,
      std::vector<Point_3>& rotated) const {

      if (angle == FT(0)) {
        
        rotated = wall;
        return;
      }

      rotated.clear();
      rotated.resize(wall.size());

      Point_3 q;
      for (std::size_t i = 0; i < wall.size(); ++i) {
        const Point_3& p = wall[i];

        q = Point_3(p.x() - b.x(), p.y() - b.y(), p.z() - b.z());
        rotate_point(angle, axis, q);
        rotated[i] = Point_3(q.x() + b.x(), q.y() + b.y(), q.z() + b.z());
      }      
    }

    void rotate_point(
      const FT angle, 
      const Vector_3& axis, 
      Point_3& p) const {

			const double tmp_angle = CGAL::to_double(angle);

			const FT c = static_cast<FT>(std::cos(tmp_angle));
			const FT s = static_cast<FT>(std::sin(tmp_angle));

			const FT C = FT(1) - c;

			const FT x = axis.x();
			const FT y = axis.y();
			const FT z = axis.z();

			p = Point_3(
        (x * x * C + c)     * p.x() + (x * y * C - z * s) * p.y() + (x * z * C + y * s) * p.z(),
				(y * x * C + z * s) * p.x() + (y * y * C + c)     * p.y() + (y * z * C - x * s) * p.z(),
				(z * x * C - y * s) * p.x() + (z * y * C + x * s) * p.y() + (z * z * C + c)     * p.z());
		}

    void triangulate_region_facets(
      const std::vector< std::vector<Point_3> >& region_facets, 
      CDT& cdt) const {

			std::vector< std::vector<Vertex_handle> > vhs;
      insert_points(region_facets, cdt, vhs);
                
      std::vector< std::pair<Vertex_handle, Vertex_handle> > final_vhs;
      update_constraints(region_facets, vhs, final_vhs);

      insert_constraints(final_vhs, cdt);
    }

    void insert_points(
      const std::vector< std::vector<Point_3> >& region_facets, 
      CDT& cdt, 
      std::vector< std::vector<Vertex_handle> >& vhs) const {
                
      cdt.clear();
      vhs.clear();

      vhs.resize(region_facets.size());
			for (std::size_t i = 0; i < region_facets.size(); ++i) {
				const auto& region_facet = region_facets[i];

				vhs[i].resize(region_facet.size());
				for (std::size_t j = 0; j < region_facet.size(); ++j) {
          const Point_3& p = region_facet[j];

					vhs[i][j] = cdt.insert(Point_2(p.x(), p.y()));
					vhs[i][j]->info().height = p.z();
				}
			}
    }

    void update_constraints(
      const std::vector< std::vector<Point_3> >& region_facets, 
      const std::vector< std::vector<Vertex_handle> >& vhs, 
      std::vector< std::pair<Vertex_handle, Vertex_handle> >& final_vhs) const {

      for (std::size_t i = 0; i < region_facets.size(); ++i) {

        for (std::size_t j = 0; j < region_facets[i].size(); ++j) {
          const std::size_t jp = (j + 1) % region_facets[i].size();

          if (is_boundary_edge(region_facets[i][j], region_facets[i][jp], i, region_facets)) {

            const auto final_constraint = std::make_pair(vhs[i][j], vhs[i][jp]);
            final_vhs.push_back(final_constraint);
          }
        }
      }
    }

    bool is_boundary_edge(
      const Point_3& p1, const Point_3& p2, 
      const std::size_t facet_index, 
      const std::vector< std::vector<Point_3> >& region_facets) const {

      for (std::size_t i = 0; i < region_facets.size(); ++i) {
        if (i == facet_index) 
          continue;

        for (std::size_t j = 0; j < region_facets[i].size(); ++j) {
          const std::size_t jp = (j + 1) % region_facets[i].size();

          if (are_equal_edges(p1, p2, region_facets[i][j], region_facets[i][jp])) 
            return false;
        }
      }
      return true;
    }

    bool are_equal_edges(
      const Point_3& p1, const Point_3& p2, 
      const Point_3& q1, const Point_3& q2) const {
      
      return 
      (are_equal_points(p1, q1) && are_equal_points(p2, q2)) || 
      (are_equal_points(p1, q2) && are_equal_points(p2, q1));
    }

    void insert_constraints(
      const std::vector< std::pair<Vertex_handle, Vertex_handle> >& final_vhs, 
      CDT& cdt) const {
                
      for (std::size_t i = 0; i < final_vhs.size(); ++i) {
        const auto& final_constraint = final_vhs[i];
                    
        if (final_constraint.first != final_constraint.second)
          cdt.insert_constraint(final_constraint.first, final_constraint.second);
      }
    }

    void get_back_region_facets(
      const CDT& cdt, 
      std::vector<Point_3>& region_facet) const {

      region_facet.clear();
      if (cdt.number_of_faces() == 0) return;

      Face_handle fh;
      bool success = find_first_face_handle(cdt, fh);
      if (!success) 
        return;

      success = traverse_cdt(fh, cdt, region_facet);
      if (!success) 
        return;

      if (region_facet.size() < 3) 
        return;
    }

    bool find_first_face_handle(
      const CDT& cdt, 
      Face_handle& fh) const {

      for (auto fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit) {
        fh = static_cast<Face_handle>(fit);

        const Vertex_handle& vh1 = fh->vertex(0);
        const Vertex_handle& vh2 = fh->vertex(1);
        const Vertex_handle& vh3 = fh->vertex(2);

        const Point_2& p1 = vh1->point();
        const Point_2& p2 = vh2->point();
        const Point_2& p3 = vh3->point();

        for (std::size_t i = 0; i < 3; ++i) {
                        
          const Edge edge = std::make_pair(fh, i);
          if (cdt.is_constrained(edge)) 
            return true;
        }
      }
      return false;
    }

    bool traverse_cdt(
      const Face_handle& fh, 
      const CDT& cdt, 
      std::vector<Point_3>& region_facet) const {
                
      Edge edge;
      region_facet.clear();

      const bool success = find_first_edge(fh, cdt, edge);
      if (!success) 
        return false;

      CGAL_precondition(edge.second >= 0 && edge.second <= 2);

      Vertex_handle vh = edge.first->vertex((edge.second + 2) % 3);
      Vertex_handle end = vh;
                
      if (vh->info().height == m_default_height) 
        return false;

      const Point_2& p = vh->point();
      region_facet.push_back(Point_3(p.x(), p.y(), vh->info().height));

      std::size_t num_iters = 0; 
      do {
        
        get_next_vertex_handle(cdt, vh, edge);
        const Point_2& q = vh->point();

        if (vh->info().height == m_default_height) 
          return false;
        
        if (vh == end) 
          break;

        region_facet.push_back(Point_3(q.x(), q.y(), vh->info().height));
                    
        if (num_iters == m_max_num_iters) 
          return false;
        ++num_iters;

      } while (vh != end);

      return is_valid_traversal(region_facet);
    }

    bool find_first_edge(
      const Face_handle& fh, 
      const CDT& cdt, 
      Edge& edge) const {

      for (int i = 0; i < 3; ++i) {              
        edge = std::make_pair(fh, i);

        if (cdt.is_constrained(edge)) 
          return true;
      }
      return false;
    }

    void get_next_vertex_handle(
      const CDT& cdt, 
      Vertex_handle& vh, 
      Edge& edge) const {

      const int index = edge.first->index(vh);
      Edge next = std::make_pair(edge.first, (index + 2) % 3);

      while (!cdt.is_constrained(next)) {

        const Face_handle fh = next.first->neighbor(next.second);
        const Vertex_handle tmp = next.first->vertex((next.second + 1) % 3);
                    
        const std::size_t tmp_index = fh->index(tmp);
        next = std::make_pair(fh, (tmp_index + 2) % 3);
      }

      vh = next.first->vertex((next.second + 2) % 3);
      edge = next;
    }

    bool is_valid_traversal(
      const std::vector<Point_3>& region_facet) const {

      if (region_facet.size() < 3) 
        return false;

      for (std::size_t i = 0; i < region_facet.size(); ++i) {
        const Point_3& p = region_facet[i];

        for (std::size_t j = 0; j < region_facet.size(); ++j) {
          const Point_3& q = region_facet[j];

          if (i == j) 
            continue;
          
          if (are_equal_points(p, q)) 
            return false;
        }
      }
      return true;
    }

    void fix_orientation(
      std::vector<Point_3>& region_facet) const {

      std::vector<Point_2> polygon(region_facet.size());
      for (std::size_t i = 0; i < region_facet.size(); ++i) {
                        
        const Point_3& p = region_facet[i];
        polygon[i] = Point_2(p.x(), p.y());
      }

      if (CGAL::orientation_2(
        polygon.begin(), polygon.end()) == CGAL::CLOCKWISE) 
        std::reverse(region_facet.begin(), region_facet.end());
    } */
  };

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_INTERNAL_WALLS_ESTIMATOR_H
