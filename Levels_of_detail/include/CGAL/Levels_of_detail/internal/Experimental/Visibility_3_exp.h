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

#ifndef CGAL_LEVELS_OF_DETAIL_INTERNAL_VISIBILITY_3_EXP_H
#define CGAL_LEVELS_OF_DETAIL_INTERNAL_VISIBILITY_3_EXP_H

// STL includes.
#include <vector>
#include <utility>

// CGAL includes.
#include <CGAL/assertions.h>

// Internal includes.
#include <CGAL/Levels_of_detail/enum.h>
#include <CGAL/Levels_of_detail/internal/struct.h>
#include <CGAL/Levels_of_detail/internal/utils.h>

// Testing.
#include "../../../../../test/Levels_of_detail/include/Saver.h"

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

  template<
  typename GeomTraits,
  typename InputRange,
  typename PointMap>
  class Visibility_3_exp {
			
  public:
    using Traits = GeomTraits;
    using Input_range = InputRange;
    using Point_map = PointMap;

    using FT = typename Traits::FT;
    using Point_2 = typename Traits::Point_2;
    using Point_3 = typename Traits::Point_3;
    using Vector_3 = typename Traits::Vector_3;
    using Line_3 = typename Traits::Line_3;
    using Plane_3 = typename Traits::Plane_3;
    using Triangle_2 = typename Traits::Triangle_2;

    using Indices = std::vector<std::size_t>;
    using Partition_3 = internal::Partition_3<Traits>;
    using Stats = std::pair<FT, FT>;
    using Face = typename Partition_3::Face;
    using Building = internal::Building<Traits>;
    
    using Triangulation = internal::Triangulation<Traits>;
    using Vertex_handle = typename Triangulation::Delaunay::Vertex_handle;

    using Saver = Saver<Traits>;

    Visibility_3_exp(
      const Input_range& input_range,
      const Point_map& point_map,
      const Building& building,
      const std::vector<Indices>& roof_points_3,
      const FT visibility_scale_3) :
    m_input_range(input_range),
    m_point_map(point_map),
    m_building(building),
    m_roof_points_3(roof_points_3),
    m_distance_tolerance(visibility_scale_3),
    m_angle_threshold(FT(10)),
    m_height_offset(m_distance_tolerance / FT(20))
    { }

    void compute(Partition_3& partition) const {
      
      if (partition.empty()) return;
      label_exterior_faces(partition.faces);

      Triangulation tri;
      preprocess(partition.faces, tri);

      for (auto& face : partition.faces) {
        if (face.exterior) {
          face.visibility = Visibility_label::OUTSIDE;
          face.inside = FT(0);
          face.outside = FT(1);
        } else compute_face_label(face);
      }
    }
    
  private:
    const Input_range& m_input_range;
    const Point_map& m_point_map;
    const Building& m_building;
    const std::vector<Indices>& m_roof_points_3;
    const FT m_distance_tolerance;
    
    // Internal parameters.
    const FT m_angle_threshold;
    const FT m_height_offset;

    void label_exterior_faces(
      std::vector<Face>& faces) const {
      
      Point_3 b;
      for (auto& face : faces) {
        face.exterior = false;

        bool exterior = false;
        const auto& neighbors = face.neighbors;
        for (const int idx : neighbors) {
          if (idx < 0) {
            exterior = true; face.exterior = exterior; break;
          }
        }
        if (exterior) continue;

        internal::compute_barycenter_3(face.vertices, b);
        if (is_above_building(b)) { 
          exterior = true; face.exterior = exterior; continue;
        }
        if (is_below_building(b)) {
          exterior = true; face.exterior = exterior; continue;
        }
        if (is_out_of_building(b)) {
          exterior = true; face.exterior = exterior; continue;
        }
      }
    }

    void compute_face_label(Face& face) const {

      face.visibility = Visibility_label::INSIDE;
      face.inside = FT(1);
      face.outside = FT(0);
    }

    bool is_above_building(const Point_3& query) const {
      return query.z() > m_building.top_z;
    }

    bool is_below_building(const Point_3& query) const {
      return query.z() < m_building.bottom_z;
    }

    bool is_out_of_building(const Point_3& query) const {
                
      const Point_2 p = Point_2(query.x(), query.y());
      const auto& tri = m_building.base1.triangulation.delaunay;
      for (auto fh = tri.finite_faces_begin();
      fh != tri.finite_faces_end(); ++fh) {
        if (!fh->info().interior) continue;
        const Triangle_2 triangle = Triangle_2(
          fh->vertex(0)->point(),
          fh->vertex(1)->point(),
          fh->vertex(2)->point());

        if (internal::is_within_triangle_2(p, triangle, m_distance_tolerance)) 
          return false;
      }
      return true;
    }

    void preprocess(
      const std::vector<Face>& polyhedrons,
      Triangulation& tri) const {

      Vector_3 v1, v2;
      std::vector<Point_3> polygon;
      tri.delaunay.clear();

      std::vector< std::vector<Point_3> > polygons;

      for (std::size_t i = 0; i < polyhedrons.size(); ++i) {
        const auto& polyhedron = polyhedrons[i];
        if (polyhedron.exterior) continue;

        const auto& faces = polyhedron.faces;
        const auto& vertices = polyhedron.vertices;
        for (const auto& face : faces) {
          
          polygon.clear();
          for (const std::size_t idx : face)
            polygon.push_back(vertices[idx]);
          const bool success = internal::compute_cross_product_3(polygon, v1);
          if (!success) std::cerr << "Error: cross product!" << std::endl;

          v2 = Vector_3(FT(0), FT(0), FT(1));
          FT angle = angle_3d(v1, v2);
          if (angle > FT(90)) angle = FT(180) - angle;
          angle = FT(90) - angle;

          if (angle < m_angle_threshold)
            continue;

          for (std::size_t k = 0; k < polygon.size(); ++k)
            tri.delaunay.insert(internal::point_2_from_point_3(polygon[k]));
        }
      }

      polygons.clear();
      polygons.resize(tri.delaunay.number_of_faces());
      std::size_t idx = 0;
      for (auto fh = tri.delaunay.finite_faces_begin();
      fh != tri.delaunay.finite_faces_end(); ++fh) {
        for (std::size_t k = 0; k < 3; ++k) {
          const Point_2& p = fh->vertex(k)->point();
          polygons[idx].push_back(Point_3(p.x(), p.y(), FT(0)));
        }
        ++idx;
      }

      Saver saver;
      saver.export_polygon_soup(polygons, "/Users/monet/Documents/lod/logs/buildings/delaunay");
    }
  };

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_INTERNAL_VISIBILITY_3_EXP_H
