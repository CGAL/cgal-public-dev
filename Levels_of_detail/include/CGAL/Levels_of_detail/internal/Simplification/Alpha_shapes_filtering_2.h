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

#define CGAL_DO_NOT_USE_BOYKOV_KOLMOGOROV_MAXFLOW_SOFTWARE
#include <CGAL/internal/Surface_mesh_segmentation/Alpha_expansion_graph_cut.h>

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
#include <CGAL/Levels_of_detail/internal/Shape_detection/Connected_component_region.h>
#include <CGAL/Levels_of_detail/internal/Shape_detection/Region_growing.h>

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
    using Vertex_handle = typename Alpha_shape_2::Vertex_handle;
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
    using Size_pair = std::pair<std::size_t, std::size_t>;
    using Alpha_expansion = CGAL::internal::Alpha_expansion_graph_cut_boost;

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

    Alpha_shapes_filtering_2(
      const FT alpha,
      const FT noise_level = FT(0)) :
    m_alpha(alpha),
    m_noise_level(noise_level),
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
        input,
        region_growing_scale_3, region_growing_angle_3,
        points);
      insert_in_triangulation(points);
    }

    void add_points_graphcut(
      const FT noise_level,
      const FT region_growing_scale_3,
      const FT region_growing_angle_3,
      const Points_3& input) {

      std::vector<Point_with_info> points;
      points.reserve(input.size());
      for (std::size_t i = 0; i < input.size(); ++i)
        points.push_back(Point_with_info(input[i], i));

      identify_wall_points(
        noise_level, region_growing_scale_3, region_growing_angle_3,
        input, points);
      insert_in_triangulation(points);
    }

    void identify_wall_points(
      const FT noise_level,
      const FT region_growing_scale_3,
      const FT region_growing_angle_3,
      const Points_3& input,
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

      // Create wall and roof points.
      Indices wall_points, roof_points;
      create_wall_and_roof_points(
        region_growing_angle_3, input, normals, wall_points, roof_points);

      // Apply region growing.
      apply_region_growing_2(
        noise_level, input, normals, wall_points);

      // Set wall indices.
      for (const std::size_t idx : wall_points)
        points[idx].belongs_to_wall = true;
    }

    void apply_region_growing_2(
      const FT noise_level,
      const Points_3& input,
      const std::vector<Vector_3>& normals,
      Indices& wall_points) {

      using Pair_item_2 = std::pair<Point_2, Vector_2>;
      using Pair_range_2 = std::vector<Pair_item_2>;
      using First_of_pair_map = CGAL::First_of_pair_property_map<Pair_item_2>;
      using Second_of_pair_map = CGAL::Second_of_pair_property_map<Pair_item_2>;

      using KNQ =
      internal::K_neighbor_query<Traits, Pair_range_2, First_of_pair_map>;
      using SNQ =
      internal::Sphere_neighbor_query<Traits, Pair_range_2, First_of_pair_map>;
      using CCR =
      internal::Connected_component_region;
      using Region_growing_2 =
      internal::Region_growing<Indices, SNQ, CCR>;

      Pair_range_2 range;
      range.reserve(wall_points.size());
      for (const std::size_t idx : wall_points) {
        const auto& p = input[idx];
        const auto& n = normals[idx];
        range.push_back(
          std::make_pair(Point_2(p.x(), p.y()), Vector_2(n.x(), n.y())));
      }

      First_of_pair_map pmap;
      SNQ neighbor_query(range, noise_level / FT(2), pmap);
      CCR region(5);

      Region_growing_2 region_growing(
        wall_points,
        neighbor_query,
        region);

      std::vector<Indices> regions;
      region_growing.detect(std::back_inserter(regions));

      Indices tmp;
      for (const auto& region : regions)
        for (const std::size_t idx : region)
          tmp.push_back(wall_points[idx]);
      wall_points = tmp;

      std::vector<Point_3> points;
      points.reserve(wall_points.size());
      for (const std::size_t idx : wall_points) {
        const auto& p = input[idx];
        points.push_back(Point_3(p.x(), p.y(), FT(0)));
      }

      Saver<Traits> saver;
      saver.export_points(
        points,
        Color(0, 0, 0),
        "/Users/monet/Documents/gf/lod/logs/buildings/tmp/wall-points-rg");
    }

    void create_wall_and_roof_points(
      const FT region_growing_angle_3,
      const Points_3& input,
      const std::vector<Vector_3>& normals,
      Indices& wall_points,
      Indices& roof_points) {

      std::vector<Point_3> save_wall_points, save_roof_points;
      const Vector_3 ref = Vector_3(FT(0), FT(0), FT(1));
      for (std::size_t i = 0; i < input.size(); ++i) {
        const auto& p = input[i];

        const auto& vec = normals[i];
        FT angle = angle_3d(vec, ref);
        if (angle > FT(90)) angle = FT(180) - angle;
        angle = FT(90) - angle;
        if (angle <= region_growing_angle_3) {
          wall_points.push_back(i);
          save_wall_points.push_back(Point_3(p.x(), p.y(), FT(0)));
        } else {
          roof_points.push_back(i);
          save_roof_points.push_back(Point_3(p.x(), p.y(), FT(0)));
        }
      }

      Saver<Traits> saver;
      saver.export_points(
        save_wall_points,
        Color(0, 0, 0),
        "/Users/monet/Documents/gf/lod/logs/buildings/tmp/wall-points");
      saver.export_points(
        save_roof_points,
        Color(0, 0, 0),
        "/Users/monet/Documents/gf/lod/logs/buildings/tmp/roof-points");
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
        if (pixel.is_interior) continue;

        const Point_2 p = Point_2(pixel.point.x(), pixel.point.y());
        Location_type type; int stub;
        const auto fh = alpha_shape.locate(p, type, stub);
        if (alpha_shape.classify(fh) == Alpha_shape_2::INTERIOR)
          pixel.is_interior = true;
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
        if (pixel.is_interior) continue;

        const Point_2 p = Point_2(pixel.point.x(), pixel.point.y());
        Location_type type; int stub;
        const auto fh = alpha_shape.locate(p, type, stub);
        if (fh->info().tagged)
          pixel.is_interior = true;
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
        "/Users/monet/Documents/gf/lod/logs/buildings/tmp/alpha_shape-original", false);

      filter_out_wrong_faces(noise_level, alpha_shape);

      save_alpha_shape(alpha_shape,
        "/Users/monet/Documents/gf/lod/logs/buildings/tmp/alpha_shape-clean", false);

      /*
      std::vector<Pair> pairs;
      for (auto fh = alpha_shape.finite_faces_begin();
      fh != alpha_shape.finite_faces_end(); ++fh) {
        if (!fh->info().tagged) continue;

        for (std::size_t k = 0; k < 3; ++k) {
          const auto& p = fh->vertex(k)->point();
          pairs.push_back(std::make_pair(p, FT(0)));
        }
      }
      m_triangulation.clear();
      insert_in_triangulation(pairs);
      set_interior_labels_stable(point_cloud); */

      for (auto& pixel : point_cloud) {
        if (pixel.is_interior) continue;

        const Point_2 p = Point_2(pixel.point.x(), pixel.point.y());
        Location_type type; int stub;
        const auto fh = alpha_shape.locate(p, type, stub);
        if (fh->info().tagged)
          pixel.is_interior = true;
      }
    }

    template<typename Pixel>
    void set_interior_labels_graphcut(
      const FT noise_level,
      std::vector<Pixel>& point_cloud) {

      if (m_triangulation.number_of_faces() == 0) return;
      CGAL_precondition(m_alpha > FT(0));
      CGAL_precondition(m_triangulation.number_of_faces() != 0);

      Alpha_shape_2 alpha_shape(
        m_triangulation, m_alpha, Alpha_shape_2::GENERAL);

      tag_faces(alpha_shape);
      save_alpha_shape(alpha_shape,
        "/Users/monet/Documents/gf/lod/logs/buildings/tmp/alpha_shape-original", false);

      use_graphcut(noise_level, alpha_shape);
      save_alpha_shape(alpha_shape,
        "/Users/monet/Documents/gf/lod/logs/buildings/tmp/alpha_shape-clean", false);

      for (auto& pixel : point_cloud) {
        if (pixel.is_interior) continue;

        const Point_2 p = Point_2(pixel.point.x(), pixel.point.y());
        Location_type type; int stub;
        const auto fh = alpha_shape.locate(p, type, stub);
        if (fh->info().tagged)
          pixel.is_interior = true;
      }
    }

  private:
    const FT m_alpha;
    const FT m_noise_level;
    Triangulation_2 m_triangulation;
    Random m_random;

    void use_graphcut(
      const FT noise_level,
      Alpha_shape_2& alpha_shape) {

      set_object_indices(alpha_shape);

      clear_labels(alpha_shape);

      label_boundary_faces(alpha_shape);

      label_wall_faces(alpha_shape);

      compute_probabilities(3, alpha_shape);

      save_alpha_shape(alpha_shape,
        "/Users/monet/Documents/gf/lod/logs/buildings/tmp/alpha_shape-init", true);

      apply_graph_cut(
        FT(1) / FT(16), true, std::size_t(-1), 3, alpha_shape);

      update_tags(alpha_shape);

      save_alpha_shape(alpha_shape,
        "/Users/monet/Documents/gf/lod/logs/buildings/tmp/alpha_shape-0", true);

      label_boundary_faces(alpha_shape);

      save_alpha_shape(alpha_shape,
        "/Users/monet/Documents/gf/lod/logs/buildings/tmp/alpha_shape-1", true);

      resize_probabilities(2, alpha_shape);

      compute_in_out(1, alpha_shape); // 1 - green - unknown
      compute_in_out(2, alpha_shape); // 2 - blue - walls
      compute_in_out(0, alpha_shape); // 0 - violet - boundary

      for (auto fh = alpha_shape.finite_faces_begin();
      fh != alpha_shape.finite_faces_end(); ++fh) {
        if (!fh->info().tagged) continue;
        const auto& probabilities = fh->info().probabilities;
        if (probabilities[1] >= FT(1) / FT(2)) // inside
          fh->info().label = 1;
        else
          fh->info().label = 0;
      }

      save_alpha_shape(alpha_shape,
        "/Users/monet/Documents/gf/lod/logs/buildings/tmp/alpha_shape-2", true);

      apply_graph_cut(
        FT(1) / FT(4), true, 0, 2, alpha_shape);

      save_alpha_shape(alpha_shape,
        "/Users/monet/Documents/gf/lod/logs/buildings/tmp/alpha_shape-3", true);

      update_tags(alpha_shape);
    }

    void set_object_indices(
      Alpha_shape_2& alpha_shape) {

      std::size_t count = 0;
      for (auto fh = alpha_shape.finite_faces_begin();
      fh != alpha_shape.finite_faces_end(); ++fh) {
        if (fh->info().tagged) {
          fh->info().object_index = count;
          ++count;
        }
      }
    }

    void clear_labels(
      Alpha_shape_2& alpha_shape) {

      for (auto fh = alpha_shape.finite_faces_begin();
      fh != alpha_shape.finite_faces_end(); ++fh) {

        if (fh->info().tagged) fh->info().label = 1;
        else fh->info().label = 0;
        fh->info().tagged_new = fh->info().tagged;
      }
    }

    void label_wall_faces(
      Alpha_shape_2& alpha_shape) {

      for (auto vh = alpha_shape.finite_vertices_begin();
      vh != alpha_shape.finite_vertices_end(); ++vh)
        if (vh->info().belongs_to_wall)
          set_incident_labels(alpha_shape, vh, 2);

      /*
      for (auto fh = alpha_shape.finite_faces_begin();
      fh != alpha_shape.finite_faces_end(); ++fh)
        if (fh->info().num_wall_points > 0)
          fh->info().label = 2; */

      /*
      using Pair_2 =
        std::pair<Point_2, std::size_t>;

      std::map<std::size_t, Face_handle> fmap;
      std::vector<Pair_2> wall_points;

      for (auto fh = alpha_shape.finite_faces_begin();
      fh != alpha_shape.finite_faces_end(); ++fh) {
        if (fh->info().label == 2) {
          const std::size_t idx = fh->info().object_index;

          fmap[idx] = fh;
          const Point_2 p = CGAL::barycenter(
            fh->vertex(0)->point(), FT(1),
            fh->vertex(1)->point(), FT(1),
            fh->vertex(2)->point(), FT(1));
          wall_points.push_back(std::make_pair(p, idx));
        }
      }

      using Pmap =
        CGAL::First_of_pair_property_map<Pair_2>;
      using SNQ =
        internal::Sphere_neighbor_query<Traits, std::vector<Pair_2>, Pmap>;

      Pmap pmap;
      SNQ snq(wall_points, m_noise_level, pmap);

      Indices neighbors;
      std::vector<Point_3> debug;
      for (std::size_t i = 0; i < wall_points.size(); ++i) {
        const auto& pair = wall_points[i];
        const auto& p = pair.first;
        snq(p, neighbors);

        if (neighbors.size() == 0)
          continue;

        const auto dense = find_dense_pair(
          wall_points, neighbors, m_noise_level / FT(2), debug);

        const auto& q = dense.first;
        const FT distance = internal::distance(p, q);

        if (distance >= m_noise_level)
          fmap[pair.second]->info().label = 1;
      }

      Saver<Traits> saver;
      saver.export_points(
        debug,
        Color(0, 0, 0),
        "/Users/monet/Documents/gf/lod/logs/buildings/tmp/dense-regions");
      */

      FT bar_area = FT(0); FT count = FT(0);
      FT max_area = FT(-1);
      for (auto fh = alpha_shape.finite_faces_begin();
      fh != alpha_shape.finite_faces_end(); ++fh) {
        if (fh->info().label != 2) continue;

        const Triangle_2 triangle = Triangle_2(
          fh->vertex(0)->point(),
          fh->vertex(1)->point(),
          fh->vertex(2)->point());
        const FT area = CGAL::abs(triangle.area());

        bar_area += area;
        count += FT(1);
        max_area = CGAL::max(area, max_area);
      }
      bar_area /= count;

      for (auto fh = alpha_shape.finite_faces_begin();
      fh != alpha_shape.finite_faces_end(); ++fh) {
        if (fh->info().label != 2) continue;

        const Triangle_2 triangle = Triangle_2(
          fh->vertex(0)->point(),
          fh->vertex(1)->point(),
          fh->vertex(2)->point());
        const FT area = CGAL::abs(triangle.area());

        const FT threshold = bar_area + (max_area - bar_area) / FT(8);
        if (area > threshold) {
          fh->info().label = 1; continue;
        }

        std::size_t num = 0;
        for (std::size_t k = 0; k < 3; ++k)
          if (fh->neighbor(k)->info().label == 0)
            num++;

        if (num != 0) {
          fh->info().label = 1; continue;
        }
      }
    }

    template<typename Pair_2>
    Pair_2 find_dense_pair(
      const std::vector<Pair_2>& wall_points,
      const Indices& neighbors,
      const FT threshold,
      std::vector<Point_3>& debug) {

      using Pmap =
        CGAL::First_of_pair_property_map<Pair_2>;
      using SNQ =
        internal::Sphere_neighbor_query<Traits, std::vector<Pair_2>, Pmap>;

      Pmap pmap;
      SNQ snq(wall_points, threshold, pmap);

      std::vector<Indices> tmp(neighbors.size());
      for (std::size_t i = 0; i < neighbors.size(); ++i) {
        const auto& pair = wall_points[neighbors[i]];
        const auto& p = pair.first;
        snq(p, tmp[i]);
      }

      int max_size = -1;
      std::size_t region_idx = std::size_t(-1);

      for (std::size_t i = 0; i < tmp.size(); ++i) {
        if (max_size < int(tmp[i].size())) {
          max_size = int(tmp[i].size()); region_idx = i;
        }
      }

      for (std::size_t idx : tmp[region_idx]) {
        const auto& p = wall_points[idx].first;
        debug.push_back(Point_3(p.x(), p.y(), FT(0)));
      }

      return wall_points[neighbors[region_idx]];
    }

    void label_boundary_faces(
      Alpha_shape_2& alpha_shape) {

      for (auto fh = alpha_shape.finite_faces_begin();
      fh != alpha_shape.finite_faces_end(); ++fh) {
        if (!fh->info().tagged) continue;

        for (std::size_t k = 0; k < 3; ++k) {
          const auto fhn = fh->neighbor(k);
          if (!fhn->info().tagged) {

            const std::size_t idx1 = (k + 1) % 3;
            const std::size_t idx2 = (k + 2) % 3;

            const auto vh1 = fh->vertex(idx1);
            const auto vh2 = fh->vertex(idx2);

            set_incident_labels(alpha_shape, vh1, 0);
            set_incident_labels(alpha_shape, vh2, 0);

            break;
          }
        }
      }
    }

    void set_incident_labels(
      const Alpha_shape_2& alpha_shape,
      const Vertex_handle vh,
      const std::size_t label) {

      auto fc = alpha_shape.incident_faces(vh);
      if (fc.is_empty()) return;
      const auto end = fc;
      do {
        if (fc->info().tagged && fc->info().label != 0)
          fc->info().label = label;
        ++fc;
      } while (fc != end);
    }

    void resize_probabilities(
      const std::size_t num_labels,
      Alpha_shape_2& alpha_shape) {

      for (auto fh = alpha_shape.finite_faces_begin();
      fh != alpha_shape.finite_faces_end(); ++fh) {
        if (!fh->info().tagged) continue;

        fh->info().probabilities.clear();
        fh->info().probabilities.resize(num_labels, FT(0));
      }
    }

    void compute_probabilities(
      const std::size_t num_labels,
      Alpha_shape_2& alpha_shape) {

      resize_probabilities(num_labels, alpha_shape);
      for (auto fh = alpha_shape.finite_faces_begin();
      fh != alpha_shape.finite_faces_end(); ++fh) {
        if (!fh->info().tagged) continue;

        for (std::size_t k = 0; k < 3; ++k) {
          const auto fhn = fh->neighbor(k);
          const std::size_t idx = fhn->info().label;
          fh->info().probabilities[idx] += FT(1);
        }

        for (std::size_t k = 0; k < num_labels; ++k)
          fh->info().probabilities[k] /= static_cast<FT>(num_labels);
      }
    }

    void update_tags(
      Alpha_shape_2& alpha_shape) {

      for (auto fh = alpha_shape.finite_faces_begin();
      fh != alpha_shape.finite_faces_end(); ++fh)
        fh->info().tagged = fh->info().tagged_new;
    }

    void compute_in_out(
      const std::size_t ref_label,
      Alpha_shape_2& alpha_shape) {

      for (auto fh = alpha_shape.finite_faces_begin();
      fh != alpha_shape.finite_faces_end(); ++fh) {
        if (!fh->info().tagged || fh->info().label != ref_label)
          continue;
        compute_statistics(
          alpha_shape, static_cast<Face_handle>(fh));
      }
    }

    void compute_statistics(
      const Alpha_shape_2& alpha_shape,
      Face_handle fh) {

      const FT radius = FT(1);
      const std::size_t num_samples = 50; // lines = num_samples / 2
      const Point_2 center = CGAL::barycenter(
        fh->vertex(0)->point(), FT(1),
        fh->vertex(1)->point(), FT(1),
        fh->vertex(2)->point(), FT(1));

      std::vector<Point_2> samples1;
      create_points_on_circle(
        center, radius, FT(0), num_samples, samples1);

      std::vector<Point_2> samples2;
      create_points_on_circle(
        center, radius, FT(180), num_samples, samples2);

      FT inside = FT(0), outside = FT(0);
      for (std::size_t k = 0; k < num_samples / 2; ++k) {
        const auto& p1 = samples1[k];
        const auto& p2 = samples2[k];

        const auto pair = get_in_out_value(
          alpha_shape, center, p1, p2, fh);

        inside  += pair.first;
        outside += pair.second;
      }

      const FT sum = inside + outside;
      inside /= sum; outside /= sum;

      auto& probabilities = fh->info().probabilities;

      probabilities[0] = outside;
      probabilities[1] = inside;
    }

    void create_points_on_circle(
      const Point_2& center,
      const FT radius,
      const FT start,
      const std::size_t num_samples,
      std::vector<Point_2>& samples) {

      samples.clear();
      samples.reserve(num_samples);

      FT factor = FT(360) / static_cast<FT>(num_samples);
      factor *= static_cast<FT>(CGAL_PI); factor /= FT(180);

      FT init = start;
      init *= static_cast<FT>(CGAL_PI); init /= FT(180);

      for (std::size_t i = 0; i < num_samples / 2; ++i) {
        const double angle =
        CGAL::to_double(init) + double(i) * CGAL::to_double(factor);

        const FT cosa = static_cast<FT>(std::cos(angle));
        const FT sina = static_cast<FT>(std::sin(angle));

        const FT x = center.x() + radius * cosa;
        const FT y = center.y() + radius * sina;

        samples.push_back(Point_2(x, y));
      }
    }

    std::pair<FT, FT> get_in_out_value(
      const Alpha_shape_2& alpha_shape,
      const Point_2& p,
      const Point_2& q1,
      const Point_2& q2,
      const Face_handle ref) {

      const auto pair1 = get_closest_label(alpha_shape, p, q1, ref);
      const auto pair2 = get_closest_label(alpha_shape, p, q2, ref);

      if (!pair1.second || !pair2.second)
        return std::make_pair(FT(1), FT(1));

      const auto fh1 = pair1.first;
      const auto fh2 = pair2.first;

      const std::size_t ref_label = ref->info().label;
      switch (ref_label) {

        case 0:
        return handle_boundary_case(ref, fh1, fh2);

        case 1:
        return handle_unknown_case(ref, fh1, fh2);

        case 2:
        return handle_wall_case(ref, fh1, fh2);

        default:
        break;
      }
      return std::make_pair(FT(1), FT(1));
    }

    std::pair<FT, FT> handle_unknown_case(
      const Face_handle ref,
      const Face_handle fh1,
      const Face_handle fh2) {

      const std::size_t l1 = fh1->info().label;
      const std::size_t l2 = fh2->info().label;

      CGAL_assertion(ref->info().label == 1);
      CGAL_assertion(l1 != 1 && l2 != 1);

      const Point_2 b1 = CGAL::barycenter(
        fh1->vertex(0)->point(), FT(1),
        fh1->vertex(1)->point(), FT(1),
        fh1->vertex(2)->point(), FT(1));

      const Point_2 b2 = CGAL::barycenter(
        fh2->vertex(0)->point(), FT(1),
        fh2->vertex(1)->point(), FT(1),
        fh2->vertex(2)->point(), FT(1));

      const FT distance = internal::distance(b1, b2);

      if (l1 == l2)
        return std::make_pair(FT(1), FT(0));

      if (distance <= m_noise_level * FT(2))
        return std::make_pair(FT(0), FT(1));

      return std::make_pair(FT(1), FT(0));
    }

    std::pair<FT, FT> handle_wall_case(
      const Face_handle ref,
      const Face_handle fh1,
      const Face_handle fh2) {

      const std::size_t l1 = fh1->info().label;
      const std::size_t l2 = fh2->info().label;

      CGAL_assertion(ref->info().label == 2);
      CGAL_assertion(l1 != 2 && l2 != 2);

      const Point_2 ref_b = CGAL::barycenter(
        ref->vertex(0)->point(), FT(1),
        ref->vertex(1)->point(), FT(1),
        ref->vertex(2)->point(), FT(1));

      const Point_2 b1 = CGAL::barycenter(
        fh1->vertex(0)->point(), FT(1),
        fh1->vertex(1)->point(), FT(1),
        fh1->vertex(2)->point(), FT(1));

      const Point_2 b2 = CGAL::barycenter(
        fh2->vertex(0)->point(), FT(1),
        fh2->vertex(1)->point(), FT(1),
        fh2->vertex(2)->point(), FT(1));

      const FT distance1 = internal::distance(ref_b, b1);
      const FT distance2 = internal::distance(ref_b, b2);

      if (l1 == 0 && l2 == 0)
        return std::make_pair(FT(1), FT(0));

      if (l1 == 1 && l2 == 1)
        return case_two_equal_labels(fh1, fh2, distance1, distance2);

      if (l1 == 1 && l2 == 0)
        return case_two_different_labels(fh1, distance1, distance2);

      if (l1 == 0 && l2 == 1)
        return case_two_different_labels(fh2, distance2, distance1);

      return std::make_pair(FT(1), FT(1));
    }

    std::pair<FT, FT> case_two_equal_labels(
      const Face_handle fh1, const Face_handle fh2,
      const FT distance1, const FT distance2) {

      const auto& probabilities1 = fh1->info().probabilities;
      const auto& probabilities2 = fh2->info().probabilities;

      const bool inside1 = ( probabilities1[1] >= FT(1) / FT(2) );
      const bool inside2 = ( probabilities2[1] >= FT(1) / FT(2) );

      if (inside1 && inside2)
        return std::make_pair(FT(1), FT(0));
      if (!inside1 && !inside2)
        return std::make_pair(FT(0), FT(1));

      if (distance1 <= distance2 && inside1)
        return std::make_pair(FT(1), FT(0));
      if (distance1  > distance2 && inside1)
        return std::make_pair(FT(0), FT(1));

      if (distance2 <= distance1 && inside2)
        return std::make_pair(FT(1), FT(0));
      if (distance2  > distance1 && inside2)
        return std::make_pair(FT(0), FT(1));

      return std::make_pair(FT(0), FT(1));
    }

    std::pair<FT, FT> case_two_different_labels(
      const Face_handle fh,
      const FT distance1, const FT distance2) {

      const auto& probabilities = fh->info().probabilities;
      const bool inside = ( probabilities[1] >= FT(1) / FT(2) );

      if (distance1 <= distance2 && inside)
        return std::make_pair(FT(1), FT(0));

      if (distance1  > distance2 && inside)
        return std::make_pair(FT(0), FT(1));

      return std::make_pair(FT(0), FT(1));
    }

    std::pair<FT, FT> case_two_different_labels(
      const Face_handle fh) {

      const auto& probabilities = fh->info().probabilities;
      const bool inside = ( probabilities[1] >= FT(1) / FT(2) );

      if (inside)
        return std::make_pair(FT(1), FT(0));

      return std::make_pair(FT(0), FT(1));
    }

    std::pair<FT, FT> handle_boundary_case(
      const Face_handle ref,
      const Face_handle fh1,
      const Face_handle fh2) {

      return std::make_pair(FT(0), FT(1)); // what about this?

      std::size_t l1 = fh1->info().label;
      std::size_t l2 = fh2->info().label;

      if (!fh1->info().tagged) l1 = 3;
      if (!fh2->info().tagged) l2 = 3;

      CGAL_assertion(ref->info().label == 0);
      CGAL_assertion(l1 != 0 && l2 != 0);

      const Point_2 ref_b = CGAL::barycenter(
        ref->vertex(0)->point(), FT(1),
        ref->vertex(1)->point(), FT(1),
        ref->vertex(2)->point(), FT(1));

      const Point_2 b1 = CGAL::barycenter(
        fh1->vertex(0)->point(), FT(1),
        fh1->vertex(1)->point(), FT(1),
        fh1->vertex(2)->point(), FT(1));

      const Point_2 b2 = CGAL::barycenter(
        fh2->vertex(0)->point(), FT(1),
        fh2->vertex(1)->point(), FT(1),
        fh2->vertex(2)->point(), FT(1));

      const FT distance1 = internal::distance(ref_b, b1);
      const FT distance2 = internal::distance(ref_b, b2);

      if (l1 == 3 && l2 == 3)
        return std::make_pair(FT(0), FT(1));

      if (l1 == 2 && l2 == 2)
        return case_two_equal_labels(fh1, fh2, distance1, distance2);

      if (l1 == 1 && l2 == 1)
        return case_two_equal_labels(fh1, fh2, distance1, distance2);

      if ( (l1 == 1 || l1 == 2) && l2 == 3)
        return case_two_different_labels(fh1);

      if (l1 == 3 && (l2 == 1 || l2 == 2))
        return case_two_different_labels(fh2);

      return std::make_pair(FT(1), FT(1));
    }

    std::pair<Face_handle, bool> get_closest_label(
      const Alpha_shape_2& alpha_shape,
      const Point_2& p,
      const Point_2& q,
      const Face_handle ref) {

      const std::size_t ref_label = ref->info().label;
      const auto error = std::make_pair(Face_handle(), false);

      LF_circulator circ = alpha_shape.line_walk(p, q, ref);
      const LF_circulator start = circ;
      if (circ.is_empty()) return error;

      do {

        if (alpha_shape.is_infinite(circ))
          return error;

        if (!circ->info().tagged && ref_label != 0)
          return error;

        if (!circ->info().tagged && ref_label == 0)
          return std::make_pair(
            static_cast<Face_handle>(circ), true);

        if (circ->info().label != ref_label)
          return std::make_pair(
            static_cast<Face_handle>(circ), true);

        ++circ;
      } while (circ != start);
      return error;
    }

    void apply_graph_cut(
      const FT beta,
      const bool use_max,
      const std::size_t ref_label,
      const std::size_t num_labels,
      Alpha_shape_2& alpha_shape) {

      std::vector<std::size_t> labels;
      set_initial_labels(alpha_shape, labels);

      std::vector<Size_pair> edges;
      std::vector<double> edge_weights;
      set_graphcut_edges(beta, use_max, alpha_shape, edges, edge_weights);

      std::vector< std::vector<double> > cost_matrix;
      set_cost_matrix(use_max, num_labels, alpha_shape, cost_matrix);

      Alpha_expansion graphcut;
      graphcut(edges, edge_weights, cost_matrix, labels);

      update_labels(ref_label, labels, alpha_shape);
    }

    void set_initial_labels(
      const Alpha_shape_2& alpha_shape,
      std::vector<std::size_t>& labels) {

      labels.clear();
      for (auto fh = alpha_shape.finite_faces_begin();
      fh != alpha_shape.finite_faces_end(); ++fh) {
        if (!fh->info().tagged) continue;
        labels.push_back(fh->info().label);
      }
    }

    void set_graphcut_edges(
      const FT beta,
      const bool use_max,
      const Alpha_shape_2& alpha_shape,
      std::vector<Size_pair>& edges,
      std::vector<double>& edge_weights) {

      edges.clear();
      edge_weights.clear();

      FT max_distance = -FT(1);
      FT sum_distance =  FT(0);

      for (auto eh = alpha_shape.finite_edges_begin();
      eh != alpha_shape.finite_edges_end(); ++eh) {
        const auto  fh = eh->first;
        const auto idx = eh->second;
        const auto fhn = fh->neighbor(idx);

        if (fh->info().tagged && fhn->info().tagged) {

          const std::size_t idxi =  fh->info().object_index;
          const std::size_t idxj = fhn->info().object_index;

          const auto& p1 = fh->vertex((idx + 1) % 3)->point();
          const auto& p2 = fh->vertex((idx + 2) % 3)->point();

          const FT distance = internal::distance(p1, p2);
          max_distance = CGAL::max(distance, max_distance);
          sum_distance += distance;

          edges.push_back(std::make_pair(idxi, idxj));
          edge_weights.push_back(CGAL::to_double(distance));
        }
      }

      if (use_max) {
        for (auto& edge_weight : edge_weights)
          edge_weight /= CGAL::to_double(max_distance);
      } else {
        for (auto& edge_weight : edge_weights)
          edge_weight /= CGAL::to_double(sum_distance);
      }

      for (auto& edge_weight : edge_weights)
        edge_weight *= CGAL::to_double(beta);
    }

    void set_cost_matrix(
      const bool use_max,
      const std::size_t num_labels,
      const Alpha_shape_2& alpha_shape,
      std::vector< std::vector<double> >& cost_matrix) {

      cost_matrix.clear();
      cost_matrix.resize(num_labels);

      FT max_area = -FT(1);
      FT sum_area =  FT(0);
      std::vector<FT> weights;

      for (auto fh = alpha_shape.finite_faces_begin();
      fh != alpha_shape.finite_faces_end(); ++fh) {
        if (!fh->info().tagged) continue;

        const auto& p0 = fh->vertex(0)->point();
        const auto& p1 = fh->vertex(1)->point();
        const auto& p2 = fh->vertex(2)->point();

        const Triangle_2 triangle = Triangle_2(p0, p1, p2);
        const FT area = CGAL::abs(triangle.area());
        max_area = CGAL::max(area, max_area);
        sum_area += area;
        weights.push_back(area);
      }

      if (use_max) {
        for (auto& weight : weights)
          weight /= max_area;
      } else {
        for (auto& weight : weights)
          weight /= sum_area;
      }

      std::size_t count = 0;
      for (auto fh = alpha_shape.finite_faces_begin();
      fh != alpha_shape.finite_faces_end(); ++fh) {
        if (!fh->info().tagged) continue;
        for (std::size_t k = 0; k < num_labels; ++k)
          cost_matrix[k].push_back(get_cost(
            weights[count], fh->info().probabilities[k]));
        ++count;
      }
    }

    double get_cost(
      const FT weight,
      const FT probability) {

      return CGAL::to_double((FT(1) - probability) * weight);
    }

    void update_labels(
      const std::size_t ref_label,
      const std::vector<std::size_t>& labels,
      Alpha_shape_2& alpha_shape) {

      std::size_t count = 0;
      for (auto fh = alpha_shape.finite_faces_begin();
      fh != alpha_shape.finite_faces_end(); ++fh) {
        if (!fh->info().tagged) continue;

        fh->info().label = labels[count];
        if (fh->info().label == ref_label)
          fh->info().tagged_new = false;
        else
          fh->info().tagged_new = true;
        ++count;
      }
    }

    void save_alpha_shape(
      const Alpha_shape_2& alpha_shape,
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
        alpha_shape, indexer, num_vertices,
        output_vertices, output_faces, z, out_labels);

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

      save_alpha_shape(alpha_shape,
        "/Users/monet/Documents/gf/lod/logs/buildings/tmp/alpha_shape-graphcut-1", true);

      retag_using_barycenter(
        noise_level, alpha_shape);

      save_alpha_shape(alpha_shape,
        "/Users/monet/Documents/gf/lod/logs/buildings/tmp/alpha_shape-tagged", true);

      set_object_indices(alpha_shape);

      compute_probabilities(3, alpha_shape);

      apply_graph_cut(
        FT(1) / FT(4), true, 0, 3, alpha_shape);

      save_alpha_shape(alpha_shape,
        "/Users/monet/Documents/gf/lod/logs/buildings/tmp/alpha_shape-graphcut-2", true);

      update_tags(alpha_shape);
    }

    void retag_using_barycenter(
      const FT noise_level,
      Alpha_shape_2& alpha_shape) {

      Point_2 b;
      compute_barycenter(alpha_shape, b);

      for (auto fh = alpha_shape.finite_faces_begin();
      fh != alpha_shape.finite_faces_end(); ++fh)
        fh->info().tagged_new = fh->info().tagged;

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

        if (found) {
          retag_along_line_using_barycenter(
            noise_level, b, alpha_shape, fh);

          for (std::size_t k = 0; k < 3; ++k) {
            const auto fhn = fh->neighbor(k);
            if (!alpha_shape.is_infinite(fhn))
              retag_along_line_using_barycenter(
                noise_level, b, alpha_shape, fhn);
          }
        }
      }

      for (auto fh = alpha_shape.finite_faces_begin();
      fh != alpha_shape.finite_faces_end(); ++fh)
        if (!fh->info().tagged_new)
          fh->info().label = 0;
    }

    void compute_barycenter(
      const Alpha_shape_2& alpha_shape,
      Point_2& b) {

      FT x = FT(0), y = FT(0), count = FT(0);
      for (auto fh = alpha_shape.finite_faces_begin();
      fh != alpha_shape.finite_faces_end(); ++fh) {
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
      const Point_2& barycenter,
      const Alpha_shape_2& alpha_shape,
      const Face_handle& fh) {

      const Point_2& p0 = fh->vertex(0)->point();
      const Point_2& p1 = fh->vertex(1)->point();
      const Point_2& p2 = fh->vertex(2)->point();

      const Point_2 b = CGAL::barycenter(
        p0, FT(1), p1, FT(1), p2, FT(1));

      const Vector_2 direction = Vector_2(FT(1), FT(0));

      const Triangle_2 triangle = Triangle_2(p0, p1, p2);
      std::vector<Point_2> samples;
      using Point_generator = CGAL::Random_points_in_triangle_2<Point_2>;

      Point_generator generator(triangle, m_random);
      std::copy_n(
        generator, 48, std::back_inserter(samples));

      for (const auto& p : samples) {
        const Line_2 line = Line_2(p, direction);
        const Point_2 q = line.projection(barycenter);
        apply_line_walk(noise_level, p, q, b, fh, alpha_shape);
      }
    }

    void apply_line_walk(
      const FT noise_level,
      const Point_2& p,
      const Point_2& q,
      const Point_2& b,
      const Face_handle& fh,
      const Alpha_shape_2& alpha_shape) {

      LF_circulator circ  = alpha_shape.line_walk(p, q, fh);
      LF_circulator start = circ;
      if (circ.is_empty()) return;

      bool found = false;
      std::size_t count = 0, num_found = 0;
      do {
        if (alpha_shape.is_infinite(circ)) break;
        if (is_closest_criteria(noise_level, b, circ)) {
          ++num_found;

          LF_circulator next = circ; ++next;
          if (!is_closest_criteria(noise_level, b, next))
            found = true;
        }
        ++circ; ++count;
      } while (circ != start && !found);

      if (count == 1)
        start->info().tagged_new = false;

      if (count > 1 && found) {
        circ = start;

        std::size_t half = 0;
        if (num_found > 1) half = std::ceil(num_found / 2);

        for (std::size_t i = 0; i < count - half; ++i) {
          circ->info().tagged_new = false;
          ++circ;
        }
      }
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
        /* if (pi.belongs_to_wall) continue; */

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

      /*
      for (std::size_t i = 0; i < points.size(); ++i) {
        const auto& pi = points[i];
        if (pi.belongs_to_wall) {

          Location_type type; int stub;
          auto fh = m_triangulation.locate(pi.point, type, stub);
          fh->info().num_wall_points += 1;
        }
      } */
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
