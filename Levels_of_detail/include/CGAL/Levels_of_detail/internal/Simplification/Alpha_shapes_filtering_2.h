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
#include <algorithm>

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

// Shape detection.
#include <CGAL/Levels_of_detail/internal/Shape_detection/Estimate_normals_3.h>

// Testing.
#include "../../../../../test/Levels_of_detail/include/Saver.h"
#include "../../../../../test/Levels_of_detail/include/Utilities.h"

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
    using Vector_2 = typename Traits::Vector_2;
    using Vector_3 = typename Traits::Vector_3;
    using Triangle_2 = typename Traits::Triangle_2;
    using Line_2 = typename Traits::Line_2;

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

    using Indices = std::vector<std::size_t>;
    using Points_2 = std::vector<Point_2>;
    using Points_3 = std::vector<Point_3>;

    using LF_circulator = typename Triangulation_2::Line_face_circulator;

    struct Point_with_info {

      Point_with_info(
        const Point_3& _point, 
        const std::size_t _idx) :
      point(Point_2(_point.x(), _point.y())),
      z(_point.z()), idx(_idx) { }
      
      const Point_2 point;
      const FT z;
      const std::size_t idx;
      
      bool belongs_to_wall = false;
    };

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

    void add_points_line_sweep(
      const FT region_growing_scale_3,
      const FT region_growing_angle_3,
      const Points_3& input) {
      
      std::vector<Point_with_info> points;
      points.reserve(input.size());
      for (std::size_t i = 0; i < input.size(); ++i)
        points.push_back(Point_with_info(input[i], i));

      identify_wall_points(
        input, region_growing_scale_3, region_growing_angle_3, 
        points);
      insert_in_triangulation(points);
    }

    void identify_wall_points(
      const Points_3& input,
      const FT region_growing_scale_3,
      const FT region_growing_angle_3,
      std::vector<Point_with_info>& points) {

      using Identity_map_3 = CGAL::Identity_property_map<Point_3>;
      using SNQ =
        internal::Sphere_neighbor_query<Traits, Points_3, Identity_map_3>;
      using NE3 = 
        internal::Estimate_normals_3<Traits, Points_3, Identity_map_3, SNQ>;

      // Compute normals.
      Identity_map_3 identity_map_3;
      std::vector<Vector_3> normals;
      SNQ neighbor_query(
        input, region_growing_scale_3, identity_map_3);
      NE3 estimator(
        input, neighbor_query, identity_map_3);
      estimator.get_normals(normals);
      CGAL_assertion(normals.size() == input.size());

      // Remove vertical points.
      std::vector<Point_3> wall_points, roof_points;
      const Vector_3 ref = Vector_3(FT(0), FT(0), FT(1));
      for (std::size_t i = 0; i < input.size(); ++i) {
        
        const auto& vec = normals[i];
        FT angle = angle_3d(vec, ref);
        if (angle > FT(90)) angle = FT(180) - angle;
        angle = FT(90) - angle;
        if (angle <= region_growing_angle_3) {
          
          points[i].belongs_to_wall = true;
          wall_points.push_back(input[i]);

        } else {
          roof_points.push_back(input[i]);
        }
      }

      Saver<Traits> saver;
      saver.export_points(
        wall_points, 
        Color(0, 0, 0),
        "/Users/monet/Documents/lod/logs/buildings/tmp/wall-points");
      saver.export_points(
        roof_points, 
        Color(0, 0, 0),
        "/Users/monet/Documents/lod/logs/buildings/tmp/roof-points");
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

      std::map<std::size_t, bool> clean;
      update_clean_points_v1(
        max_height_difference, bounds, indices, pairs, clean);

      /*
      update_clean_points_v2(
        noise_level, bounds, indices, pairs, clean); */

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

      Pmap pmap;
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

      Alpha_shape_2 alpha_shape(
        m_triangulation, m_alpha, Alpha_shape_2::GENERAL);

      for (auto& pixel : point_cloud) {
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

      Alpha_shape_2 alpha_shape(
        m_triangulation, m_alpha, Alpha_shape_2::GENERAL);
      tag_faces(alpha_shape);

      for (auto& pixel : point_cloud) {
        const Point_2 p = Point_2(pixel.point.x(), pixel.point.y());
        Location_type type; int stub;
        const auto fh = alpha_shape.locate(p, type, stub);
        if (fh->info().tagged)
          pixel.is_interior = true;
        else 
          pixel.is_interior = false;
      }
    }

    template<typename Pixel>
    void set_interior_labels_line_sweep(
      const FT noise_level,
      std::vector<Pixel>& point_cloud) {

      if (m_triangulation.number_of_faces() == 0) return;
      CGAL_precondition(m_alpha > FT(0));
      CGAL_precondition(m_triangulation.number_of_faces() != 0);

      Alpha_shape_2 alpha_shape(
        m_triangulation, m_alpha, Alpha_shape_2::GENERAL);

      tag_faces(alpha_shape);

      save_alpha_shape(alpha_shape, 
        "/Users/monet/Documents/lod/logs/buildings/tmp/alpha_shape-original", false);

      filter_out_wrong_faces(noise_level, alpha_shape);

      save_alpha_shape(alpha_shape, 
        "/Users/monet/Documents/lod/logs/buildings/tmp/alpha_shape-tagged", true);

      for (auto& pixel : point_cloud) {
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

    void save_alpha_shape(
      const Alpha_shape_2& tri,
      const std::string path,
      const bool out_labels) {

      const FT z = FT(0);
      std::size_t num_vertices = 0;
      internal::Indexer<Point_3> indexer;

      std::vector<Point_3> vertices; 
      std::vector<Indices> faces; 
      std::vector<Color> fcolors;

      Polygon_inserter<Traits> inserter(faces, fcolors);
      auto output_vertices = std::back_inserter(vertices);
      auto output_faces = boost::make_function_output_iterator(inserter);

      output_triangulation(
        tri, indexer, num_vertices, output_vertices, output_faces, z, out_labels);
      
      Saver<Traits> saver;
      saver.export_polygon_soup(vertices, faces, fcolors, path);
    }

    template<
    typename Indexer,
    typename VerticesOutputIterator,
    typename FacesOutputIterator>
    void output_triangulation(
      const Alpha_shape_2& tri,
      Indexer& indexer,
      std::size_t& num_vertices,
      VerticesOutputIterator vertices,
      FacesOutputIterator faces,
      const FT z,
      const bool out_labels) const {

      std::vector<std::size_t> face(3);
      for (auto fh = tri.finite_faces_begin(); 
      fh != tri.finite_faces_end(); ++fh) {
        if (!fh->info().tagged) continue;
        
        for (std::size_t k = 0; k < 3; ++k) {
          const Point_2& q = fh->vertex(k)->point();
          const Point_3 p = Point_3(q.x(), q.y(), z);
          const std::size_t idx = indexer(p);
          if (idx == num_vertices) {
            *(vertices++) = p; 
            ++num_vertices;
          }
          face[k] = idx;
        }
        if (out_labels)
          *(faces++) = std::make_pair(face, fh->info().label);
        else 
          *(faces++) = std::make_pair(face, 1);
      }
    }

    void filter_out_wrong_faces(
      const FT noise_level,
      Alpha_shape_2& alpha_shape) {

      retag_using_barycenter(noise_level, alpha_shape);
      for (auto fh = alpha_shape.finite_faces_begin();
      fh != alpha_shape.finite_faces_end(); ++fh)
        fh->info().tagged = fh->info().tagged_new;
    }

    void retag_using_barycenter(
      const FT noise_level,
      Alpha_shape_2& alpha_shape) {

      Point_2 b;
      compute_barycenter(alpha_shape, b);

      for (auto fh = alpha_shape.finite_faces_begin();
      fh != alpha_shape.finite_faces_end(); ++fh) {
        if (!fh->info().tagged) continue;

        bool found = false;
        for (std::size_t k = 0; k < 3; ++k) {
          const auto fhn = fh->neighbor(k);
          if (!fhn->info().tagged) {
            found = true; break;
          }
        }

        if (found)
          retag_along_line_using_barycenter(
            noise_level, b, alpha_shape, fh);
      }
    }

    void compute_barycenter(
      const Alpha_shape_2& alpha_shape,
      Point_2& b) {
      
      FT x = FT(0), y = FT(0), count = FT(0);
      for (auto fh = alpha_shape.finite_faces_begin();
      fh != alpha_shape.finite_faces_end(); ++fh) {
        fh->info().tagged_new = fh->info().tagged;
        if (!fh->info().tagged) continue;

        const Point_2 p = CGAL::barycenter(
          fh->vertex(0)->point(), FT(1),
          fh->vertex(1)->point(), FT(1),
          fh->vertex(2)->point(), FT(1));

        x += p.x();
        y += p.y();

        count += FT(1);
      }
      x /= count;
      y /= count;

      b = Point_2(x, y);
    }

    void retag_along_line_using_barycenter(
      const FT noise_level,
      const Point_2& b,
      const Alpha_shape_2& alpha_shape, 
      const Face_handle& fh) {
      
      const Point_2& p0 = fh->vertex(0)->point();
      const Point_2& p1 = fh->vertex(1)->point();
      const Point_2& p2 = fh->vertex(2)->point();
      
      const Vector_2 direction = Vector_2(FT(1), FT(0));

      const Triangle_2 triangle = Triangle_2(p0, p1, p2);
      std::vector<Point_2> samples;
      using Point_generator = CGAL::Random_points_in_triangle_2<Point_2>;

      Point_generator generator(triangle, m_random);
      std::copy_n(
        generator, 12, std::back_inserter(samples));

      for (const auto& p : samples) {
        const Line_2 line = Line_2(p, direction);
        const Point_2 q = line.projection(b);
        apply_line_walk(noise_level, p, q, fh, alpha_shape);
      }
    }

    void apply_line_walk(
      const FT noise_level,
      const Point_2& p,
      const Point_2& q,
      const Face_handle& fh,
      const Alpha_shape_2& alpha_shape) {

      LF_circulator circ  = alpha_shape.line_walk(p, q, fh);
      LF_circulator start = circ;
      if (circ.is_empty()) return;

      LF_circulator stop;
      bool found = false; std::size_t count = 0;
      do {
        if (alpha_shape.is_infinite(circ)) break;
        if (is_closest_criteria(noise_level, p, circ)) {
          stop = circ; found = true;
        }
        ++circ; ++count; 
      } while (circ != start && !found);

      if (count == 1) {
        start->info().tagged_new = false;
      }

      if (count > 1 && found) {
        circ = start; 
        circ->info().tagged_new = false;
        do {
          ++circ; 
          circ->info().tagged_new = false;
        } while (circ != stop);
      }

      fh->info().label = 3; 
    }

    bool is_closest_criteria(
      const FT noise_level,
      const Point_2& p,
      LF_circulator& circ) {
      
      if (circ->info().label == 2) {
          
        const Point_2 q = CGAL::barycenter(
          circ->vertex(0)->point(), FT(1),
          circ->vertex(1)->point(), FT(1),
          circ->vertex(2)->point(), FT(1));

        const FT dist = internal::distance(p, q);
        if (dist <= noise_level * FT(2))
          return true;
      }
      return false;
    }

    void update_clean_points_v2(
      const FT noise_level,
      const std::vector<Pair>& bounds,
      const Indices& indices,
      std::vector<Pair>& pairs,
      std::map<std::size_t, bool>& clean) {

      Pmap pmap;
      Sphere_neighbor_query neighbor_query(
        pairs, noise_level, pmap);

      for (std::size_t i = 0; i < pairs.size(); ++i)
        clean[i] = false;

      Indices neighbors;
      for (std::size_t i = 0; i < bounds.size(); ++i) {
        const auto& bound = bounds[i];
        
        const auto& query = bound.first;
        neighbor_query(query, neighbors);
        neighbors.push_back(indices[i]);
        handle_neighborhood(noise_level, neighbors, pairs, clean);
      }
    }

    void handle_neighborhood(
      const FT noise_level,
      const Indices& neighbors,
      std::vector<Pair>& pairs,
      std::map<std::size_t, bool>& clean) {

      for (const std::size_t idx : neighbors)
        clean[idx] = true;
      return;

      Indices local_indices;
      std::vector<Pair> local_pairs;

      local_pairs.reserve(neighbors.size());
      local_indices.reserve(neighbors.size());

      for (const std::size_t idx : neighbors) {
        local_pairs.push_back(pairs[idx]);
        local_indices.push_back(idx);
      }

      Pmap pmap;
      Sphere_neighbor_query neighbor_query(
        local_pairs, noise_level / FT(4), pmap);

      std::vector<Indices> local_neighbors;
      local_neighbors.reserve(local_pairs.size());

      for (const auto& local_pair : local_pairs) {
        const auto& query = local_pair.first;

        Indices tmp;
        neighbor_query(query, tmp);
        local_neighbors.push_back(tmp);
      }

      std::sort(local_neighbors.begin(), local_neighbors.end(), 
      [](const Indices& a, const Indices& b) -> bool { 
        return a.size() > b.size();
      });

      Indices final_indices;
      final_indices.push_back(0);

      for (const std::size_t idx : neighbors)
        clean[idx] = false;
      for (const std::size_t final_index : final_indices)
        for (const std::size_t idx : local_neighbors[final_index])
          clean[local_indices[idx]] = true;
    }

    void update_clean_points_v1(
      const FT max_height_difference,
      const std::vector<Pair>& bounds,
      const Indices& indices,
      const std::vector<Pair>& pairs,
      std::map<std::size_t, bool>& clean) {

      Pmap pmap;
      K_neighbor_query neighbor_query(
        pairs, FT(6), pmap);

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

    void tag_faces(
      Alpha_shape_2& alpha_shape) {

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

    void insert_in_triangulation(
      const std::vector<Point_with_info>& points) {

      m_triangulation.clear();
      for (std::size_t i = 0; i < points.size(); ++i) {
        const auto& pi = points[i];

        auto vh = m_triangulation.insert(pi.point);
        vh->info().z = pi.z;
        vh->info().object_index = pi.idx;
        vh->info().belongs_to_wall = pi.belongs_to_wall;
      }

      for (auto fh = m_triangulation.finite_faces_begin();
      fh != m_triangulation.finite_faces_end(); ++fh)
        fh->info().label = 1;

      for (auto vh = m_triangulation.finite_vertices_begin();
      vh != m_triangulation.finite_vertices_end(); ++vh) {
        if (vh->info().belongs_to_wall) {
          
          auto fc = m_triangulation.incident_faces(vh);
          if (fc.is_empty()) continue;
          const auto end = fc;
          do {
            fc->info().label = 2; ++fc;
          } while (fc != end);
        }
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
