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

#ifndef CGAL_LEVELS_OF_DETAIL_INTERNAL_VISIBILITY_3_H
#define CGAL_LEVELS_OF_DETAIL_INTERNAL_VISIBILITY_3_H

// STL includes.
#include <vector>
#include <utility>

// CGAL includes.
#include <CGAL/assertions.h>
#include <CGAL/number_utils.h>
#include <CGAL/Polygon_2_algorithms.h>
#include <CGAL/Barycentric_coordinates_2/Mean_value_2.h>
#include <CGAL/Barycentric_coordinates_2/Segment_coordinates_2.h>
#include <CGAL/Barycentric_coordinates_2/Generalized_barycentric_coordinates_2.h>
#include <CGAL/Delaunay_triangulation_3.h>

// Internal includes.
#include <CGAL/Levels_of_detail/enum.h>
#include <CGAL/Levels_of_detail/internal/struct.h>

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

  template<
  typename GeomTraits,
  typename InputRange,
  typename PointMap>
  class Visibility_3 {
			
  public:
    using Traits = GeomTraits;
    using Input_range = InputRange;
    using Point_map = PointMap;

    using FT = typename Traits::FT;

    using Indices = std::vector<std::size_t>;
    using Partition_3 = internal::Partition_3<Traits>;

    /*
    using Point_2 = typename Traits::Point_2;
    using Point_3 = typename Traits::Point_3;
    using Triangle_2 = typename Traits::Triangle_2;
    using Line_2 = typename Traits::Line_2;
    using Line_3 = typename Traits::Line_3;
    using Vector_3 = typename Traits::Vector_3;
    using Plane_3 = typename Traits::Plane_3;
    
    using Building = Building<Traits>;
    using Polyhedron_facet_3 = Polyhedron_facet_3<Traits>;
    
    using Stats = std::pair<FT, FT>;

    typename Traits::Compute_squared_distance_2 squared_distance_2;
    typename Traits::Construct_cross_product_vector_3 cross_product_3;
    typename Traits::Compute_squared_length_3 squared_length_3;
    typename Traits::Compute_scalar_product_3 dot_product_3;

    using Mean_value = 
      CGAL::Barycentric_coordinates::Mean_value_2<Traits>;
    using Mean_value_coordinates = 
      CGAL::Barycentric_coordinates::Generalized_barycentric_coordinates_2<Mean_value, Traits>;

    using Intersect = typename Traits::Intersect_3;

    using Delaunay_3 = CGAL::Delaunay_triangulation_3<Traits>; */

    Visibility_3(
      const Input_range& input_range,
      const Point_map& point_map,
      const std::vector<Indices>& roof_points_3) :
    m_input_range(input_range),
    m_point_map(point_map),
    m_roof_points_3(roof_points_3)
    /*,
    m_building_max_height(FT(0)),
    m_building_min_height(FT(0)),
    m_distance_tolerance(FT(2)),
    m_bc_tolerance_top(FT(6) / FT(5)),
    m_bc_tolerance_bottom(-FT(1) / FT(5)),
    m_tolerance(FT(1) / FT(100000)),
    m_angle_threshold(FT(10)),
    m_height_offset(FT(1) / FT(8)) */ { 

      // compute_building_max_height();
      // compute_building_min_height();
    }

    void compute(Partition_3& partition) const {
      
      /*
      for (std::size_t i = 0; i < polyhedrons.size(); ++i) {
        auto& polyhedron = polyhedrons[i];

        const Stats stats = estimate_in_out_values(polyhedron);

        CGAL_precondition(stats.first  >= FT(0) && stats.first  <= FT(1));
        CGAL_precondition(stats.second >= FT(0) && stats.second <= FT(1));
                    
        CGAL_precondition(
          CGAL::abs(stats.first + stats.second - FT(1)) < FT(1) / FT(100000));

        polyhedron.in = stats.first;
        polyhedron.out = stats.second;

        polyhedron.weight = compute_weight(polyhedron);

        if (polyhedron.in > FT(1) / FT(2)) 
          polyhedron.visibility = Visibility_label::INSIDE;
        else 
          polyhedron.visibility = Visibility_label::OUTSIDE;
      }
      normalize_weights(polyhedrons); */
    }
    
  private:
    const Input_range& m_input_range;
    const Point_map& m_point_map;
    const std::vector<Indices>& m_roof_points_3;
    
    /*
    FT m_building_max_height;
    FT m_building_min_height;

    const FT m_distance_tolerance;
    const FT m_bc_tolerance_top;
    const FT m_bc_tolerance_bottom;

    const FT m_tolerance;
    const FT m_angle_threshold;
    const FT m_height_offset; */

    /*
    void compute_building_max_height() {

      // Alternatively I can use the whole m_input_range here
      // as in the original code!

      const auto& roof_indices = m_building.roof_indices;
      CGAL_precondition(roof_indices.size() > 0);

      m_building_max_height = -internal::max_value<FT>();
			for (std::size_t i = 0; i < roof_indices.size(); ++i) {
        
        const auto& indices = roof_indices[i];
        for (std::size_t j = 0; j < indices.size(); ++j) {
          
          const Point_3& p = 
            get(m_point_map, *(m_input_range.begin() + indices[j]));
          
          m_building_max_height = 
            CGAL::max(m_building_max_height, p.z());
        }
      }
    }

    void compute_building_min_height() {

      const auto& ground = m_building.approximate_ground;
      CGAL_precondition(ground.size() > 0);

      m_building_min_height = internal::max_value<FT>();
			for (std::size_t i = 0; i < ground.size(); ++i) {
        const Point_3& p = ground[i];
        
        m_building_min_height = 
          CGAL::min(m_building_min_height, p.z());
      }
    }

    Stats estimate_in_out_values(const Polyhedron_facet_3& polyhedron) const {

      Point_3 barycenter;
      compute_polyhedron_barycenter(polyhedron, barycenter);

      if (is_above_building_max_height(barycenter)) 
        return std::make_pair(FT(0), FT(1));

      if (is_below_building_min_height(barycenter)) 
        return std::make_pair(FT(0), FT(1));

      if (is_out_of_building(barycenter)) 
        return std::make_pair(FT(1) / FT(5), FT(4) / FT(5));

      if (has_vertices_outside(polyhedron)) 
        return std::make_pair(FT(2) / FT(5), FT(3) / FT(5));
                
      return estimate_in_out_values_statistically(polyhedron, barycenter);
    }

    FT compute_weight(const Polyhedron_facet_3& polyhedron) const {
      return compute_volume(polyhedron);
    }

    void normalize_weights(
      std::vector<Polyhedron_facet_3>& polyhedrons) const {

			FT total_weight = FT(0);
			for (std::size_t i = 0; i < polyhedrons.size(); ++i) {
					
				const auto& polyhedron = polyhedrons[i];
				total_weight += polyhedron.weight;
			}

			for (std::size_t i = 0; i < polyhedrons.size(); ++i)
				polyhedrons[i].weight /= total_weight;
    }

    void compute_polyhedron_barycenter(
      const Polyhedron_facet_3& polyhedron, 
      Point_3& b) const {

      const std::vector<Point_3>& vertices = polyhedron.vertices;
      const std::size_t num_vertices = vertices.size();

      CGAL_precondition(num_vertices > 0);

      FT x = FT(0), y = FT(0), z = FT(0);
      for (std::size_t i = 0; i < num_vertices; ++i) {

        x += vertices[i].x();
        y += vertices[i].y();
        z += vertices[i].z();
      }

      x /= static_cast<FT>(num_vertices);
      y /= static_cast<FT>(num_vertices);
      z /= static_cast<FT>(num_vertices);

      b = Point_3(x, y, z);
    }

    bool is_above_building_max_height(const Point_3& query) const {
      return query.z() > m_building_max_height;
    }

    bool is_below_building_min_height(const Point_3& query) const {
      return query.z() < m_building_min_height;
    }

    bool is_out_of_building(const Point_3& query) const {
                
      const Point_2 p = Point_2(query.x(), query.y());

      const auto& footprint = m_building.triangles;
      for (std::size_t i = 0; i < footprint.size(); ++i) {              
        const Triangle_2& triangle = footprint[i];

        if (is_within_triangle(p, triangle)) 
          return false;
      }
      return true;
    }

    bool is_within_triangle(
      const Point_2& query, 
      const Triangle_2& triangle) const {
                
      if (triangle.has_on_bounded_side(query) || 
          triangle.has_on_boundary(query)) 
        return true;
                
      for (std::size_t i = 0; i < 3; ++i) {
        const std::size_t ip = (i + 1) % 3;

        const Point_2& p1 = triangle.vertex(i);
        const Point_2& p2 = triangle.vertex(ip);

        const Line_2 line = Line_2(p1, p2);

        const Point_2 projected = line.projection(query);
        const FT squared_distance = squared_distance_2(query, projected);

        const CGAL::cpp11::array<FT, 2> pair = 
          Barycentric_coordinates::
          compute_segment_coordinates_2(p1, p2, projected, Traits());

        const FT squared_tolerance = 
          m_distance_tolerance * m_distance_tolerance;
                    
        const FT epst = m_bc_tolerance_top;
        const FT epsb = m_bc_tolerance_bottom;

        if (pair[0] > epsb && 
            pair[1] > epsb && 
            pair[0] < epst && 
            pair[1] < epst && 
            squared_distance < squared_tolerance) 
          return true;
      }
      return false;
    }

    bool has_vertices_outside(const Polyhedron_facet_3& polyhedron) const {

      std::size_t count = 0;
      const auto& vertices = polyhedron.vertices;
      for (std::size_t i = 0; i < vertices.size(); ++i) {

        const auto& vertex = vertices[i];
        const bool is_out = is_out_of_building(vertex);

        if (is_out) 
          ++count;
        
        if (is_out && count > 0) 
          return true;
      }
      return false;
    }

    Stats estimate_in_out_values_statistically(
      const Polyhedron_facet_3& polyhedron,
      const Point_3& barycenter) const {

      const auto& vertices = polyhedron.vertices;
      const auto& faces = polyhedron.faces;

      std::size_t in = 0, out = 0;
      std::vector<std::size_t> indices;

      for (std::size_t i = 0; i < faces.size(); ++i) {    
        const auto& face = faces[i];

        if (is_vertical_face(vertices, face)) 
          continue;
        
        process_face(vertices, face, indices, in, out);
      }
      process_middle_plane(barycenter, indices, in, out);
                
      if (in == 0 && out == 0) 
        return std::make_pair(FT(1) / FT(5), FT(4) / FT(5));

      const FT tmp_in = static_cast<FT>(in);
      const FT tmp_out = static_cast<FT>(out);

      const FT sum = tmp_in + tmp_out;
      CGAL_precondition(sum != FT(0));

      const FT final_in = tmp_in  / sum;
      const FT final_out = tmp_out / sum;

      return std::make_pair(final_in, final_out);
    }

    bool is_vertical_face(
      const std::vector<Point_3>& vertices, 
      const std::vector<std::size_t>& face) const {

			Vector_3 face_normal;
			const bool success = set_face_normal(vertices, face, face_normal);

      if (!success) 
        return true;

			const Vector_3 ground_normal = Vector_3(FT(0), FT(0), FT(1));

      const FT angle = compute_angle(face_normal, ground_normal);
      const FT angle_diff = CGAL::abs(FT(90) - CGAL::abs(angle));

      if (angle_diff < m_angle_threshold) 
        return true;
    
      return false;
    }

    bool set_face_normal(
      const std::vector<Point_3>& vertices, 
      const std::vector<std::size_t>& face, 
      Vector_3& face_normal) const {

      CGAL_precondition(face.size() >= 3);

      const Vector_3 zero = Vector_3(FT(0), FT(0), FT(0));
      for (std::size_t i = 0; i < face.size(); ++i) {

        const std::size_t ip = (i + 1) % face.size();
        const std::size_t ipp = (i + 2) % face.size();
                    
        const Point_3& p1 = vertices[face[i]];
        const Point_3& p2 = vertices[face[ip]];
        const Point_3& p3 = vertices[face[ipp]];

        const Vector_3 v1 = Vector_3(p2, p1);
        const Vector_3 v2 = Vector_3(p2, p3);

        face_normal = cross_product_3(v1, v2);
        if (!are_equal_points(face_normal, zero)) {
                     
          normalize(face_normal);
          return true;
        }
      }
      return false;
		}

    template<class Point>
    bool are_equal_points(const Point& p, const Point& q) const {

      const FT eps = m_tolerance;
      return 
        (CGAL::abs(p.x() - q.x()) < eps) && 
        (CGAL::abs(p.y() - q.y()) < eps) && 
        (CGAL::abs(p.z() - q.z()) < eps);
    }

    void normalize(Vector_3& v) const {
      v /= static_cast<FT>(
        CGAL::sqrt(
          CGAL::to_double(
            v.squared_length())));
    }

    FT compute_angle(const Vector_3& m, const Vector_3& n) const {
				
			const auto cross = cross_product_3(m, n);
			const FT length = static_cast<FT>(
        CGAL::sqrt(
          CGAL::to_double(
            squared_length_3(cross))));

			const FT dot = dot_product_3(m, n);

			FT angle_rad = static_cast<FT>(
        std::atan2(
          CGAL::to_double(length), 
          CGAL::to_double(dot)));
                
      const FT half_pi = static_cast<FT>(CGAL_PI) / FT(2);
      
      if (angle_rad > half_pi) 
        angle_rad = static_cast<FT>(CGAL_PI) - angle_rad;

			const FT angle_deg = 
        angle_rad * FT(180) / static_cast<FT>(CGAL_PI);
      
      return angle_deg;
		}

    void process_face(
      const std::vector<Point_3>& vertices, 
      const std::vector<std::size_t>& face, 
      std::vector<std::size_t>& indices, 
      std::size_t& in, 
      std::size_t& out) const {

      std::vector<Point_2> polygon;
      create_polygon(vertices, face, polygon);

      if (!CGAL::is_simple_2(polygon.begin(), polygon.end())) 
        return;

      Mean_value_coordinates mvc(polygon.begin(), polygon.end());

      // Alternatively I can use the whole m_input_range here
      // as in the original code!
      const auto& roof_indices = m_building.roof_indices;
      for (std::size_t i = 0; i < roof_indices.size(); ++i) {
        for (std::size_t j = 0; j < roof_indices[i].size(); ++j) {

          const std::size_t point_index = roof_indices[i][j];
          const Point_3& p = 
            get(m_point_map, *(m_input_range.begin() + point_index));

          const FT height = 
            intersect_face(mvc, vertices, face, p);

          if (height == internal::max_value<FT>()) 
            continue;

          indices.push_back(point_index);
          
          if (is_inside_building(height, p.z())) 
            ++in;
          else 
            ++out;
        }
      }
    }

    void create_polygon(
      const std::vector<Point_3>& vertices, 
      const std::vector<std::size_t>& face, 
      std::vector<Point_2>& polygon) const {

      polygon.clear();
      polygon.resize(face.size());

      for (std::size_t i = 0; i < face.size(); ++i) {
                    
        const Point_3& p = vertices[face[i]];
        polygon[i] = Point_2(p.x(), p.y());
      }
    }

    FT intersect_face(
      Mean_value_coordinates& mvc, 
      const std::vector<Point_3>& vertices, 
      const std::vector<std::size_t>& face, 
      const Point_3& p) const {

      const Point_2 query = Point_2(p.x(), p.y());

      std::vector<FT> coordinates;
      mvc(query, std::back_inserter(coordinates));
            
      if (is_inside_polygon(coordinates)) 
        return intersect_line_and_plane(vertices, face, p);
      
      return internal::max_value<FT>();
    }

    bool is_inside_polygon(const std::vector<FT>& coordinates) const {

      CGAL_precondition(coordinates.size() >= 3);
      for (std::size_t i = 0 ; i < coordinates.size(); ++i)
        if (coordinates[i] <= FT(0) || coordinates[i] >= FT(1)) 
          return false;
      return true;
    }

    FT intersect_line_and_plane(
      const std::vector<Point_3>& vertices, 
      const std::vector<std::size_t>& face, 
      const Point_3& p) const {

      Line_3 line;
      create_line(p, line);

      Plane_3 plane;
      const bool success = create_plane(vertices, face, plane);

      if (!success) 
        return internal::max_value<FT>();

      return intersect(line, plane);
    }

    void create_line(const Point_3& p, Line_3& line) const {
                
      const Point_3 p1 = Point_3(p.x(), p.y(), m_building_min_height - FT(10));
      const Point_3 p2 = Point_3(p.x(), p.y(), m_building_max_height + FT(10));

      line = Line_3(p1, p2);
    }

    bool create_plane(
      const std::vector<Point_3>& vertices, 
      const std::vector<std::size_t>& face, 
      Plane_3& plane) const {

      CGAL_precondition(face.size() >= 3);

      const Vector_3 zero = Vector_3(FT(0), FT(0), FT(0));
      for (std::size_t i = 0; i < face.size(); ++i) {

        const std::size_t ip = (i + 1) % face.size();
        const std::size_t ipp = (i + 2) % face.size();
                    
        const Point_3& p1 = vertices[face[i]];
        const Point_3& p2 = vertices[face[ip]];
        const Point_3& p3 = vertices[face[ipp]];

        const Vector_3 v1 = Vector_3(p2, p1);
        const Vector_3 v2 = Vector_3(p2, p3);

        const Vector_3 normal = cross_product_3(v1, v2);
        if (!are_equal_points(normal, zero)) {
                     
          plane = Plane_3(p1, p2, p3);
          return true;
        }
      }
      return false;
    }

    FT intersect(const Line_3& line, const Plane_3& plane) const {

			typename CGAL::cpp11::result_of<Intersect(Line_3, Plane_3)>::type result 
      = intersection(line, plane);
                
      if (result) {
        if (const Line_3* tmp = boost::get<Line_3>(&*result)) 
          return internal::max_value<FT>();
        else {
          const Point_3* point = boost::get<Point_3>(&*result);
				  return (*point).z();
        }
      }
      return internal::max_value<FT>();
    }

    bool is_inside_building(
      const FT current_height, 
      const FT real_height) const {
      
      return 
        current_height > m_building_min_height - m_height_offset && 
        current_height < real_height + m_height_offset;
    }

    void process_middle_plane(
      const Point_3& barycenter, 
      const std::vector<std::size_t>& indices, 
      size_t& in, 
      size_t& out) const {

      Plane_3 middle_plane(
        barycenter, Vector_3(FT(0), FT(0), FT(1)));

      Line_3 line;
      for (std::size_t i = 0; i < indices.size(); ++i) {
                    
        const std::size_t point_index = indices[i];
        const Point_3& p = 
          get(m_point_map, *(m_input_range.begin() + point_index));

        create_line(p, line);
        const FT height = intersect(line, middle_plane);

        if (is_inside_building(height, p.z())) 
          ++in;
        else 
          ++out;
      }
    }

    FT compute_volume(const Polyhedron_facet_3& polyhedron) const {
                
      Delaunay_3 delaunay_3;
      create_triangulation_3(polyhedron, delaunay_3);

      FT total_volume = FT(0);
      for (auto it = delaunay_3.finite_cells_begin(); 
      it != delaunay_3.finite_cells_end(); ++it) {
        
        const auto& tetrahedron = delaunay_3.tetrahedron(it);

        const FT volume = tetrahedron.volume();
        total_volume += volume;
      }
      return total_volume;
    }

    void create_triangulation_3(
      const Polyhedron_facet_3& polyhedron, 
      Delaunay_3& delaunay_3) const {

      delaunay_3.clear();

      const auto& vertices = polyhedron.vertices;
      for (std::size_t i = 0; i < vertices.size(); ++i) {
                    
        const Point_3& p = vertices[i];
        delaunay_3.insert(p);
      }
    } */
  };

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_INTERNAL_VISIBILITY_3_H
