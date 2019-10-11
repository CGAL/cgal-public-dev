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

#ifndef CGAL_LEVELS_OF_DETAIL_INTERNAL_ALPHA_SHAPES_FILTERING_2_H
#define CGAL_LEVELS_OF_DETAIL_INTERNAL_ALPHA_SHAPES_FILTERING_2_H

// STL includes.
#include <map>
#include <vector>

// CGAL includes.
#include <CGAL/Alpha_shape_2.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/Random.h>

// Internal includes.
#include <CGAL/Levels_of_detail/internal/utils.h>
#include <CGAL/Levels_of_detail/internal/struct.h>

// Spatial search.
#include <CGAL/Levels_of_detail/internal/Spatial_search/K_neighbor_query.h>
#include <CGAL/Levels_of_detail/internal/Spatial_search/Sphere_neighbor_query.h>

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

  template<typename GeomTraits>
  class Alpha_shapes_filtering_2 {

  public:
    using Traits = GeomTraits;
    using FT = typename Traits::FT;
    using Point_2 = typename Traits::Point_2;
    using Point_3 = typename Traits::Point_3;
    using Triangle_2 = typename Traits::Triangle_2;

    using Fi = Face_info<Traits>;
    using Fbi = CGAL::Triangulation_face_base_with_info_2<Fi, Traits>;
    using Fb = CGAL::Alpha_shape_face_base_2<Traits, Fbi>;

    using Vi = Vertex_info<Traits>;
    using Vbi = CGAL::Triangulation_vertex_base_with_info_2<Vi, Traits>;
    using Vb = CGAL::Alpha_shape_vertex_base_2<Traits, Vbi>;
    
    using Tds = CGAL::Triangulation_data_structure_2<Vb, Fb>;
    using Triangulation_2 = CGAL::Delaunay_triangulation_2<Traits, Tds>;
    using Alpha_shape_2 = CGAL::Alpha_shape_2<Triangulation_2>;
    using Points_in_triangle = CGAL::Random_points_in_triangle_2<Point_2>;
    using Location_type = typename Triangulation_2::Locate_type;
    using Face_handle = typename Alpha_shape_2::Face_handle;
    using Random = CGAL::Random;

    using Pair = std::pair<Point_2, FT>;
    using Pmap = CGAL::First_of_pair_property_map<Pair>;

    using K_neighbor_query =
      internal::K_neighbor_query<Traits, std::vector<Pair>, Pmap>;
    using Sphere_neighbor_query =
      internal::Sphere_neighbor_query<Traits, std::vector<Pair>, Pmap>;
    using Neighbor_query = 
      K_neighbor_query;

    using Indices = std::vector<std::size_t>;
    using Points_2 = std::vector<Point_2>;
    using Points_3 = std::vector<Point_3>;

    Alpha_shapes_filtering_2(const FT alpha) : 
    m_alpha(alpha),
    m_random(0) 
    { }

    template<
    typename Range, 
    typename Point_map>
    void add_points(const Range& range, Point_map point_map) {
      insert_in_triangulation(range, point_map);
    }

    void add_points_with_filtering(
      const FT noise_level,
      const FT max_height_difference,
      const Points_3& points,
      Points_3& result) {
      
      std::vector<Pair> pairs;
      pairs.reserve(points.size());
      for (const auto& p : points) 
        pairs.push_back(
          std::make_pair(Point_2(p.x(), p.y()), p.z()));
      
      insert_in_triangulation(pairs);
      Alpha_shape_2 alpha_shape(
        m_triangulation, m_alpha, Alpha_shape_2::GENERAL);
      tag_faces(alpha_shape);

      std::vector<Pair> bounds; Indices indices;
      get_bounds_of_alpha_shapes(alpha_shape, bounds, indices);

      Pmap pmap;
      Neighbor_query neighbor_query(
        pairs, FT(6), pmap);

      std::map<std::size_t, bool> clean;
      update_clean_points_v1(
        max_height_difference, bounds, indices, pairs, neighbor_query, clean);

      result.clear(); std::vector<Pair> finals;
      for (const auto& pair : clean) {
        if (pair.second) {
          const auto& p = pairs[pair.first];
          result.push_back(
            Point_3(p.first.x(), p.first.y(), p.second));
          finals.push_back(
            std::make_pair(Point_2(p.first.x(), p.first.y()), p.second));
        }
      }

      insert_in_triangulation(finals, pmap);
    }

    void get_filtered_points(
      const FT edge_sampling_2, 
      Points_2& result) {
      
      if (m_triangulation.number_of_faces() == 0) return;
      CGAL_precondition(m_alpha > FT(0));
      CGAL_precondition(m_triangulation.number_of_faces() != 0);

      Alpha_shape_2 alpha_shape(m_triangulation, m_alpha, Alpha_shape_2::GENERAL);
      sample_edges(alpha_shape, edge_sampling_2, result);
    }

    void get_samples(
      const FT edge_sampling_2,
      const std::size_t face_sampling_2, 
      Points_2& result) {
      
      if (m_triangulation.number_of_faces() == 0) return;
      CGAL_precondition(m_alpha > FT(0));
      CGAL_precondition(m_triangulation.number_of_faces() != 0);

      Alpha_shape_2 alpha_shape(m_triangulation, m_alpha, Alpha_shape_2::GENERAL);
      sample_faces(alpha_shape, edge_sampling_2, face_sampling_2, result);
    }

    template<typename Pixel>
    void set_interior_labels_stable(
      std::vector<Pixel>& point_cloud) {

      if (m_triangulation.number_of_faces() == 0) return;
      CGAL_precondition(m_alpha > FT(0));
      CGAL_precondition(m_triangulation.number_of_faces() != 0);

      Alpha_shape_2 alpha_shape(m_triangulation, m_alpha, Alpha_shape_2::GENERAL);

      for (auto& pixel : point_cloud) {
        /* if (pixel.is_interior) continue; */
        
        const Point_2 p = Point_2(pixel.point.x(), pixel.point.y());
        Location_type type; int stub;
        const auto fh = alpha_shape.locate(p, type, stub);
        if (alpha_shape.classify(fh) == Alpha_shape_2::INTERIOR)
          pixel.is_interior = true;
        else 
          pixel.is_interior = false;
      }
    }

    template<typename Pixel>
    void set_interior_labels_tagged(
      std::vector<Pixel>& point_cloud) {

      if (m_triangulation.number_of_faces() == 0) return;
      CGAL_precondition(m_alpha > FT(0));
      CGAL_precondition(m_triangulation.number_of_faces() != 0);

      Alpha_shape_2 alpha_shape(m_triangulation, m_alpha, Alpha_shape_2::GENERAL);
      tag_faces(alpha_shape);

      for (auto& pixel : point_cloud) {
        /* if (pixel.is_interior) continue; */
        
        const Point_2 p = Point_2(pixel.point.x(), pixel.point.y());
        Location_type type; int stub;
        const auto fh = alpha_shape.locate(p, type, stub);
        if (fh->info().tagged)
          pixel.is_interior = true;
        else 
          pixel.is_interior = false;
      }
    }

  private:
    const FT m_alpha;
    Triangulation_2 m_triangulation;
    Random m_random;

    void update_clean_points_v2(
      const std::vector<Pair>& bounds,
      const std::vector<Pair>& pairs,
      Neighbor_query& neighbor_query,
      std::map<std::size_t, bool>& clean) {

      for (std::size_t i = 0; i < pairs.size(); ++i)
        clean[i] = true;

      
    }

    void update_clean_points_v1(
      const FT max_height_difference,
      const std::vector<Pair>& bounds,
      const Indices& indices,
      const std::vector<Pair>& pairs,
      Neighbor_query& neighbor_query,
      std::map<std::size_t, bool>& clean) {

      for (std::size_t i = 0; i < pairs.size(); ++i)
        clean[i] = true;

      Indices neighbors;
      for (std::size_t i = 0; i < bounds.size(); ++i) {
        const auto& a = bounds[i];
        
        const auto& ap = a.first;
        neighbor_query(ap, neighbors);
        const FT az = a.second;

        FT maxz = -FT(1);
        for (const std::size_t idx : neighbors) {
          const auto& b = pairs[idx];
          const FT bz = b.second;

          maxz = CGAL::max(maxz, bz);
        }
        maxz = CGAL::max(maxz, az);

        for (const std::size_t idx : neighbors) {
          const auto& b = pairs[idx];
          const FT bz = b.second;

          const FT diff = CGAL::abs(maxz - bz);
          if (diff > max_height_difference)
            clean[idx] = false;
        }

        const FT diff = CGAL::abs(maxz - az);
        if (diff > max_height_difference)
          clean[indices[i]] = false;
      }
    }

    void get_bounds_of_alpha_shapes(
      const Alpha_shape_2& alpha_shape,
      std::vector<Pair>& bounds,
      Indices& indices) {

      bounds.clear(); indices.clear();
      for (auto fh = alpha_shape.finite_faces_begin();
      fh != alpha_shape.finite_faces_end(); ++fh) {
    
        for (std::size_t k = 0; k < 3; ++k) {
          const auto fhn = fh->neighbor(k);
          if (alpha_shape.is_infinite(fhn))
            continue;

          if ( 
            ( fh->info().tagged && !fhn->info().tagged) || 
            (fhn->info().tagged &&  !fh->info().tagged) ) {

            const auto vh = fh->vertex( (k + 1) % 3 );
            const auto& p = vh->point();

            indices.push_back(vh->info().object_index);
            bounds.push_back(std::make_pair(p, vh->info().z));
          }
        }
      }
    }

    void tag_faces(Alpha_shape_2& alpha_shape) {

      for (auto fit = alpha_shape.finite_faces_begin();
      fit != alpha_shape.finite_faces_end(); ++fit)
        fit->info().tagged = true;

      for (auto fh = alpha_shape.finite_faces_begin();
      fh != alpha_shape.finite_faces_end(); ++fh) {
        
        bool found = false;
        for (std::size_t k = 0; k < 3; ++k) {
          const auto fhn = fh->neighbor(k);
          if (alpha_shape.is_infinite(fhn)) {
            found = true; break;
          }
        }
        if (!found) continue;
        if (!fh->info().tagged) continue;

        Face_handle seed = static_cast<Face_handle>(fh);
        propagate(alpha_shape, seed);
      }
    }

    void propagate(
      const Alpha_shape_2& alpha_shape, 
      Face_handle& fh) {

      if (alpha_shape.classify(fh) == Alpha_shape_2::INTERIOR)
        return;
      fh->info().tagged = false;

      for (std::size_t k = 0; k < 3; ++k) {
        auto fhn = fh->neighbor(k);
        if (!alpha_shape.is_infinite(fhn) && fhn->info().tagged)
          propagate(alpha_shape, fhn);
      }
    }

    void insert_in_triangulation(
      const std::vector<Pair>& range) {

      m_triangulation.clear();  
      for (std::size_t i = 0; i < range.size(); ++i) {
        auto vh = m_triangulation.insert(
          internal::point_2_from_point_3(range[i].first));
        vh->info().z = range[i].second;
        vh->info().object_index = i;
      }
    }

    template<
    typename Range, 
    typename Point_map>
    void insert_in_triangulation(
      const Range& range, 
      Point_map point_map) {
                  
      for (auto it = range.begin(); it != range.end(); ++it)
        m_triangulation.insert(
          internal::point_2_from_point_3(get(point_map, *it)));
    }

    void sample_edges(
      const Alpha_shape_2& alpha_shape,
      const FT edge_sampling_2,
      Points_2& result) const {

      for (auto eit = alpha_shape.alpha_shape_edges_begin(); 
        eit != alpha_shape.alpha_shape_edges_end(); ++eit) {

        const Point_2& source = eit->first->vertex((eit->second + 1) % 3)->point();
        const Point_2& target = eit->first->vertex((eit->second + 2) % 3)->point();
        sample_edge(source, target, edge_sampling_2, result);
      }
    }

    void sample_edge(
      const Point_2& source, 
      const Point_2& target,
      const FT edge_sampling_2,
      Points_2& result) const {

      const FT distance = internal::distance(source, target);
        
      CGAL_precondition(edge_sampling_2 > FT(0));
      const std::size_t nb_pts = 
      static_cast<std::size_t>(CGAL::to_double(distance / edge_sampling_2)) + 1;
        
      CGAL_precondition(nb_pts > 0);
      for (std::size_t i = 0; i <= nb_pts; ++i) {

        const FT ratio = static_cast<FT>(i) / static_cast<FT>(nb_pts);
        result.push_back(
          Point_2(
            source.x() * (FT(1) - ratio) + target.x() * ratio,
            source.y() * (FT(1) - ratio) + target.y() * ratio));
      }
    }

    void sample_faces(
      const Alpha_shape_2& alpha_shape,
      const FT edge_sampling_2,
      const std::size_t face_sampling_2,
      Points_2& result) {

      for (auto fit = alpha_shape.finite_faces_begin(); 
        fit != alpha_shape.finite_faces_end(); ++fit) {

        const auto type = alpha_shape.classify(fit);
        if (type == Alpha_shape_2::INTERIOR) {
          
          const auto& p0 = fit->vertex(0)->point();
          const auto& p1 = fit->vertex(1)->point();
          const auto& p2 = fit->vertex(2)->point();

          // Add edges.
          sample_edge(p0, p1, edge_sampling_2, result);
          sample_edge(p1, p2, edge_sampling_2, result);
          sample_edge(p2, p0, edge_sampling_2, result);

          // Add face.
          const Triangle_2 triangle(p0, p1, p2);
          Points_in_triangle generator(triangle, m_random);
          std::copy_n(generator, face_sampling_2, 
          std::back_inserter(result));
        }
      }
    }
  };

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_INTERNAL_ALPHA_SHAPES_FILTERING_2_H
