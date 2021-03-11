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

#ifndef CGAL_LEVELS_OF_DETAIL_INTERNAL_PARTITION_BUILDER_FROM_IMAGE_2_H
#define CGAL_LEVELS_OF_DETAIL_INTERNAL_PARTITION_BUILDER_FROM_IMAGE_2_H

// STL includes.
#include <map>
#include <set>
#include <list>
#include <queue>
#include <vector>
#include <utility>
#include <memory>

// CGAL includes.
#include <CGAL/property_map.h>
#include <CGAL/Polyline_simplification_2/simplify.h>

// Internal includes.
#include <CGAL/Levels_of_detail/internal/utils.h>
#include <CGAL/Levels_of_detail/internal/struct.h>

// Other includes.
#include <CGAL/Levels_of_detail/internal/Shape_detection/Region_growing.h>
#include <CGAL/Levels_of_detail/internal/Shape_detection/Linear_image_region.h>
#include <CGAL/Levels_of_detail/internal/Shape_detection/Planar_image_region.h>
#include <CGAL/Levels_of_detail/internal/Shape_detection/Oriented_image_region.h>
#include <CGAL/Levels_of_detail/internal/Spatial_search/Oriented_neighbor_query.h>
#include <CGAL/Levels_of_detail/internal/Regularization/Polygon_regularizer.h>
#include <CGAL/Levels_of_detail/internal/Spatial_search/K_neighbor_query.h>

// Testing.
#include "../../../../../test/Levels_of_detail/include/Saver.h"
#include "../../../../../test/Levels_of_detail/include/Utilities.h"

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

template<
typename GeomTraits,
typename ImagePointer>
class Partition_builder_from_image_2 {

public:
  using Traits = GeomTraits;

  using FT          = typename Traits::FT;
  using Point_2     = typename Traits::Point_2;
  using Point_3     = typename Traits::Point_3;
  using Vector_2    = typename Traits::Vector_2;
  using Vector_3    = typename Traits::Vector_3;
  using Line_2      = typename Traits::Line_2;
  using Line_3      = typename Traits::Line_3;
  using Segment_2   = typename Traits::Segment_2;
  using Segment_3   = typename Traits::Segment_3;
  using Plane_3     = typename Traits::Plane_3;
  using Intersect_2 = typename Traits::Intersect_2;
  using Intersect_3 = typename Traits::Intersect_3;

  using Triangulation    = internal::Triangulation<Traits>;
  using Face_handle      = typename Triangulation::Delaunay::Face_handle;
  using Vertex_handle    = typename Triangulation::Delaunay::Vertex_handle;
  using LF_circulator    = typename Triangulation::Delaunay::Line_face_circulator;
  using Partition_edge_2 = internal::Partition_edge_2<Traits>;
  using Partition_face_2 = internal::Partition_face_2<Traits>;

  using Partition_2 = internal::Partition_2<Traits>;
  using Size_pair   = std::pair<std::size_t, std::size_t>;
  using Idx_map     = std::map<Size_pair, std::size_t>;
  using Indices     = std::vector<std::size_t>;
  using FT_pair     = std::pair<FT, FT>;
  using Seg_pair    = std::pair<Segment_2, bool>;

  using Saver = Saver<Traits>;

  using Point_pair = std::pair<Point_2, std::size_t>;
  using PS_generator = CGAL::Points_on_segment_2<Point_2>;
  using Point_pair_map = CGAL::First_of_pair_property_map<Point_pair>;
  using KNQ_pair =
    internal::K_neighbor_query<Traits, std::vector<Point_pair>, Point_pair_map>;
  using Polygon_regularizer = internal::Polygon_regularizer<Traits>;

  enum class Point_type {
    DEFAULT = 0,
    UNIQUE_FREE = 1,
    UNIQUE_LINEAR = 2,
    FREE = 3,
    LINEAR = 4,
    CORNER = 5,
    BOUNDARY = 6
  };

  struct Pixel {
    Point_2 point;
    Indices neighbors_03;
    Indices neighbors_47;
    bool is_outer = false;
    std::size_t index = std::size_t(-1);
    std::size_t label = std::size_t(-1);
    std::size_t i = std::size_t(-1), j = std::size_t(-1);
    std::size_t binary = std::size_t(-1);
    bool used = false;
  };

  struct My_segment {
    Point_2 source_;
    Point_2 target_;
    std::set<std::size_t> ls, lt;

    const Point_2& source() const {
      return source_;
    }
    const Point_2& target() const {
      return target_;
    }
  };

  struct Segment_wrapper {
    std::size_t index = std::size_t(-1);
    Indices neighbors;
  };

  struct Regular_segment {
    Point_2 source_;
    Point_2 target_;
    bool skip = false;
    Point_type source_type = Point_type::DEFAULT;
    Point_type target_type = Point_type::DEFAULT;
    Point_type saved_source_type = Point_type::DEFAULT;
    Point_type saved_target_type = Point_type::DEFAULT;
    std::size_t source_bd_idx = std::size_t(-1);
    std::size_t target_bd_idx = std::size_t(-1);

    const Point_2& source() const {
      return source_;
    }
    const Point_2& target() const {
      return target_;
    }
  };

  struct My_point {
    Point_2 point_;
    std::set<std::size_t> labels;
    bool belongs_to_line = false;
    Point_type type = Point_type::DEFAULT;
    Point_type saved_type = Point_type::DEFAULT;
    std::size_t bd_idx = std::size_t(-1);

    const Point_2& point() const {
      return point_;
    }
  };

  struct Contour {
    std::vector<My_point> points;
    std::vector<My_point> simplified;
    std::vector<My_point> regularized;
    std::set<Size_pair> neighbors;
    std::set<Size_pair> ns0;
    std::set<Size_pair> ns1;
    bool is_closed = false;
  };

  struct Image {
    std::vector<Pixel> pixels;
    Indices seeds;
    std::vector<Size_pair> label_pairs;
    std::size_t num_labels;
    std::vector<My_segment> segments;
    std::vector<Contour> contours;
    std::vector<Pixel> dual;
    Line_2 direction;

    bool is_ridge = false;

    FT noise_level;
    FT min_length_2;
    FT angle_bound_2;
    FT ordinate_bound_2;

    void use_version_4() {
      m_use_version_8 = false;
    }

    void use_version_8() {
      m_use_version_8 = true;
    }

    void operator()(
      const std::size_t query_index,
      Indices& neighbors) {

      neighbors.clear();
      CGAL_assertion(query_index >= 0 && query_index < pixels.size());

      const auto& ns03 = pixels[query_index].neighbors_03;
      for (const std::size_t neighbor : ns03)
        if (neighbor != std::size_t(-1))
          neighbors.push_back(neighbor);

      if (m_use_version_8) {
        const auto& ns47 = pixels[query_index].neighbors_47;
        for (const std::size_t neighbor : ns47)
          if (neighbor != std::size_t(-1))
            neighbors.push_back(neighbor);
      }
    }

    void clear() {
      pixels.clear();
      seeds.clear();
      label_pairs.clear();
      segments.clear();
      contours.clear();
      dual.clear();
      num_labels = 0;
    }

    void create_contours() {
      if (!is_ridge) return;

      make_binary_indices();
      make_dual_grid();
      apply_contouring();
      make_contours();
    }

    void simplify_contours(
      const std::map<std::size_t, Plane_3>& plane_map) {

      const bool line_found = intersect_labels(plane_map, direction);
      for (auto& contour : contours) {
        if (line_found)
          set_linear_points(direction, contour);
        simplify_contour(direction, contour);
      }
    }

    void regularize_contours(
      const std::vector<Segment_2>& boundaries,
      const std::vector<Point_pair>& boundary_queries,
      std::shared_ptr<KNQ_pair>& knq_ptr) {

      for (auto& contour : contours)
        regularize_contour(
          boundaries, boundary_queries, knq_ptr, contour);
    }

  private:
    bool m_use_version_8 = false;

    void make_binary_indices() {

      CGAL_assertion(label_pairs.size() == 1);
      const std::size_t ref_label = label_pairs[0].first;
      for (auto& pixel : pixels) {
        if (pixel.label == ref_label)
          pixel.binary = 0;
        else
          pixel.binary = 1;
      }
    }

    void make_dual_grid() {

      /* save_original_grid(pixels,
        Color(125, 0, 0), Color(0, 0, 125)); */

      std::vector<Point_2> points;
      for (auto& pixel1 : pixels) {
        const auto& neighbors03 = pixel1.neighbors_03;
        const auto& neighbors47 = pixel1.neighbors_47;

        std::size_t neighbors03_size = 0;
        for (const std::size_t neighbor : neighbors03)
          if (neighbor != std::size_t(-1))
            neighbors03_size += 1;

        std::size_t neighbors47_size = 0;
        for (const std::size_t neighbor : neighbors47)
          if (neighbor != std::size_t(-1))
            neighbors47_size += 1;

        if (
          neighbors03_size < 4 ||
          neighbors47_size < 4) continue;
        pixel1.used = true;

        CGAL_assertion(neighbors47.size() == 4);
        for (std::size_t i = 0; i < neighbors47.size(); ++i) {
          const std::size_t neighbor = neighbors47[i];
          if (neighbor == std::size_t(-1)) continue;

          const auto& pixel2 = pixels[neighbor];
          if (pixel2.used) continue;

          const auto& p = pixel1.point;
          const auto& q = pixel2.point;
          const auto  m = internal::middle_point_2(p, q);
          points.push_back(m);
        }
      }

      std::sort(points.begin(), points.end());
      points.erase(
        std::unique(
          points.begin(), points.end(), [](const Point_2& p, const Point_2& q) {
          return internal::are_equal_points_2(p, q);
        }),
      points.end());

      dual.clear();
      dual.reserve(points.size());

      Pixel pixel;
      for (const auto& point : points) {
        pixel.point = point;
        pixel.neighbors_47.clear();
        pixel.neighbors_47.resize(4, std::size_t(-1));
        dual.push_back(pixel);
      }

      for (const auto& pixel1 : pixels) {
        const auto& neighbors = pixel1.neighbors_47;
        for (std::size_t i = 0; i < neighbors.size(); ++i) {

          const std::size_t neighbor = neighbors[i];
          if (neighbor == std::size_t(-1)) continue;
          const auto& pixel2 = pixels[neighbor];

          const auto& p = pixel1.point;
          const auto& q = pixel2.point;
          const auto  m = internal::middle_point_2(p, q);

          for (auto& pixel : dual) {
            if (internal::are_equal_points_2(pixel.point, m)) {
              if (i == 0) {
                pixel.neighbors_47[0] = pixel2.index;
                pixel.neighbors_47[2] = pixel1.index;
              }

              if (i == 1) {
                pixel.neighbors_47[1] = pixel2.index;
                pixel.neighbors_47[3] = pixel1.index;
              }

              if (i == 2) {
                pixel.neighbors_47[0] = pixel1.index;
                pixel.neighbors_47[2] = pixel2.index;
              }

              if (i == 3) {
                pixel.neighbors_47[1] = pixel1.index;
                pixel.neighbors_47[3] = pixel2.index;
              }
              break;
            }
          }
        }
      }
      /* save_dual_grid(dual, Color(0, 125, 0)); */
    }

    void apply_contouring() {

      segments.clear();
      for (const auto& pixel : dual) {
        const auto& neighbors = pixel.neighbors_47;
        const std::size_t cell_idx = get_cell_idx(neighbors);
        if (cell_idx == std::size_t(-1)) continue;

        switch (cell_idx) {
          case 0:  { add_segment_case0(neighbors);  break; }
          case 1:  { add_segment_case1(neighbors);  break; }
          case 2:  { add_segment_case2(neighbors);  break; }
          case 3:  { add_segment_case3(neighbors);  break; }
          case 4:  { add_segment_case4(neighbors);  break; }
          case 5:  { add_segment_case5(neighbors);  break; }
          case 6:  { add_segment_case6(neighbors);  break; }
          case 7:  { add_segment_case7(neighbors);  break; }
          case 8:  { add_segment_case8(neighbors);  break; }
          case 9:  { add_segment_case9(neighbors);  break; }
          case 10: { add_segment_case10(neighbors); break; }
          case 11: { add_segment_case11(neighbors); break; }
          case 12: { add_segment_case12(neighbors); break; }
          case 13: { add_segment_case13(neighbors); break; }
          default : break;
        }
      }

      /*
      std::vector<Segment_2> saved;
      saved.reserve(segments.size());
      for (const auto& segment : segments)
        saved.push_back(Segment_2(segment.source(), segment.target()));

      Saver saver;
      saver.save_polylines(
        saved, "/Users/monet/Documents/gf/lod/logs/buildings/tmp/ms_contour"); */
    }

    std::size_t get_cell_idx(const Indices& ns) {
      CGAL_assertion(ns.size() == 4);

      const std::size_t i0 = ns[0];
      const std::size_t i1 = ns[1];
      const std::size_t i2 = ns[2];
      const std::size_t i3 = ns[3];

      CGAL_assertion(i0 != std::size_t(-1));
      CGAL_assertion(i1 != std::size_t(-1));
      CGAL_assertion(i2 != std::size_t(-1));
      CGAL_assertion(i3 != std::size_t(-1));

      const std::size_t b0 = pixels[i0].binary;
      const std::size_t b1 = pixels[i1].binary;
      const std::size_t b2 = pixels[i2].binary;
      const std::size_t b3 = pixels[i3].binary;

      CGAL_assertion(b0 != std::size_t(-1));
      CGAL_assertion(b1 != std::size_t(-1));
      CGAL_assertion(b2 != std::size_t(-1));
      CGAL_assertion(b3 != std::size_t(-1));

      /* std::cout << i0 << " " << i1 << " " << i2 << " " << i3 << std::endl; */

      if (b0 == 0 && b1 == 0 && b2 == 0 && b3 == 0)
        return std::size_t(-1);
      if (b0 == 1 && b1 == 1 && b2 == 1 && b3 == 1)
        return std::size_t(-1);

      if (b0 == 0 && b1 == 1 && b2 == 1 && b3 == 1)
        return 0;
      if (b0 == 1 && b1 == 0 && b2 == 1 && b3 == 1)
        return 1;
      if (b0 == 1 && b1 == 1 && b2 == 0 && b3 == 1)
        return 2;
      if (b0 == 1 && b1 == 1 && b2 == 1 && b3 == 0)
        return 3;

      if (b0 == 1 && b1 == 0 && b2 == 0 && b3 == 0)
        return 4;
      if (b0 == 0 && b1 == 1 && b2 == 0 && b3 == 0)
        return 5;
      if (b0 == 0 && b1 == 0 && b2 == 1 && b3 == 0)
        return 6;
      if (b0 == 0 && b1 == 0 && b2 == 0 && b3 == 1)
        return 7;

      if (b0 == 0 && b1 == 0 && b2 == 1 && b3 == 1)
        return 8;
      if (b0 == 1 && b1 == 0 && b2 == 0 && b3 == 1)
        return 9;
      if (b0 == 1 && b1 == 1 && b2 == 0 && b3 == 0)
        return 10;
      if (b0 == 0 && b1 == 1 && b2 == 1 && b3 == 0)
        return 11;

      if (b0 == 0 && b1 == 1 && b2 == 0 && b3 == 1)
        return 12;
      if (b0 == 1 && b1 == 0 && b2 == 1 && b3 == 0)
        return 13;

      return std::size_t(-1);
    }

    void add_segment(
      const Pixel& px0, const Pixel& px1, const Pixel& px2) {

      My_segment segment;
      segment.source_ = internal::middle_point_2(px0.point, px1.point);
      segment.target_ = internal::middle_point_2(px1.point, px2.point);
      segment.ls.insert(px0.label);
      segment.ls.insert(px1.label);
      segment.lt.insert(px1.label);
      segment.lt.insert(px2.label);
      segments.push_back(segment);
    }

    void add_segment(
      const Pixel& px0, const Pixel& px1, const Pixel& px2, const Pixel& px3) {

      My_segment segment;
      segment.source_ = internal::middle_point_2(px0.point, px1.point);
      segment.target_ = internal::middle_point_2(px2.point, px3.point);
      segment.ls.insert(px0.label);
      segment.ls.insert(px1.label);
      segment.lt.insert(px2.label);
      segment.lt.insert(px3.label);
      segments.push_back(segment);
    }

    void add_segment_case0(const Indices& ns) {
      const auto& px0 = pixels[ns[3]];
      const auto& px1 = pixels[ns[0]];
      const auto& px2 = pixels[ns[1]];
      add_segment(px0, px1, px2);
    }

    void add_segment_case1(const Indices& ns) {
      const auto& px0 = pixels[ns[0]];
      const auto& px1 = pixels[ns[1]];
      const auto& px2 = pixels[ns[2]];
      add_segment(px0, px1, px2);
    }

    void add_segment_case2(const Indices& ns) {
      const auto& px0 = pixels[ns[1]];
      const auto& px1 = pixels[ns[2]];
      const auto& px2 = pixels[ns[3]];
      add_segment(px0, px1, px2);
    }

    void add_segment_case3(const Indices& ns) {
      const auto& px0 = pixels[ns[2]];
      const auto& px1 = pixels[ns[3]];
      const auto& px2 = pixels[ns[0]];
      add_segment(px0, px1, px2);
    }

    void add_segment_case4(const Indices& ns) {
      const auto& px0 = pixels[ns[3]];
      const auto& px1 = pixels[ns[0]];
      const auto& px2 = pixels[ns[1]];
      add_segment(px0, px1, px2);
    }

    void add_segment_case5(const Indices& ns) {
      const auto& px0 = pixels[ns[0]];
      const auto& px1 = pixels[ns[1]];
      const auto& px2 = pixels[ns[2]];
      add_segment(px0, px1, px2);
    }

    void add_segment_case6(const Indices& ns) {
      const auto& px0 = pixels[ns[1]];
      const auto& px1 = pixels[ns[2]];
      const auto& px2 = pixels[ns[3]];
      add_segment(px0, px1, px2);
    }

    void add_segment_case7(const Indices& ns) {
      const auto& px0 = pixels[ns[2]];
      const auto& px1 = pixels[ns[3]];
      const auto& px2 = pixels[ns[0]];
      add_segment(px0, px1, px2);
    }

    void add_segment_case8(const Indices& ns) {
      const auto& px0 = pixels[ns[3]];
      const auto& px1 = pixels[ns[0]];
      const auto& px2 = pixels[ns[1]];
      const auto& px3 = pixels[ns[2]];
      add_segment(px0, px1, px2, px3);
    }

    void add_segment_case9(const Indices& ns) {
      const auto& px0 = pixels[ns[2]];
      const auto& px1 = pixels[ns[3]];
      const auto& px2 = pixels[ns[0]];
      const auto& px3 = pixels[ns[1]];
      add_segment(px0, px1, px2, px3);
    }

    void add_segment_case10(const Indices& ns) {
      const auto& px0 = pixels[ns[3]];
      const auto& px1 = pixels[ns[0]];
      const auto& px2 = pixels[ns[1]];
      const auto& px3 = pixels[ns[2]];
      add_segment(px0, px1, px2, px3);
    }

    void add_segment_case11(const Indices& ns) {
      const auto& px0 = pixels[ns[2]];
      const auto& px1 = pixels[ns[3]];
      const auto& px2 = pixels[ns[0]];
      const auto& px3 = pixels[ns[1]];
      add_segment(px0, px1, px2, px3);
    }

    void add_segment_case12(const Indices& ns) {
      const auto& px0 = pixels[ns[0]];
      const auto& px1 = pixels[ns[1]];
      const auto& px2 = pixels[ns[2]];
      add_segment(px0, px1, px2);
      const auto& px3 = pixels[ns[2]];
      const auto& px4 = pixels[ns[3]];
      const auto& px5 = pixels[ns[0]];
      add_segment(px3, px4, px5);
    }

    void add_segment_case13(const Indices& ns) {
      const auto& px0 = pixels[ns[3]];
      const auto& px1 = pixels[ns[0]];
      const auto& px2 = pixels[ns[1]];
      add_segment(px0, px1, px2);
      const auto& px3 = pixels[ns[1]];
      const auto& px4 = pixels[ns[2]];
      const auto& px5 = pixels[ns[3]];
      add_segment(px3, px4, px5);
    }

    void make_contours() {
      contours.clear();
      std::vector<Segment_wrapper> segs;
      create_segment_wrappers(segs);

      using ONQ = internal::Oriented_neighbor_query<Traits, Segment_wrapper>;
      using OIR = internal::Oriented_image_region<Traits, Segment_wrapper>;
      using Region_growing = internal::Region_growing<
        std::vector<Segment_wrapper>, ONQ, OIR, Seed_map>;

      Indices seeds, idx_map;
      seeds.resize(segs.size());
      idx_map.resize(segs.size());

      for (std::size_t i = 0; i < segs.size(); ++i) {
        seeds[i] = segs[i].index; idx_map[segs[i].index] = i;
      }

      Seed_map seed_map(seeds);
      ONQ onq(segs, idx_map);
      OIR oir(segs, idx_map);

      Linear_image_region linear_region;
      Region_growing region_growing(
        segs, onq, oir, seed_map);

      std::vector<Indices> regions;
      region_growing.detect(std::back_inserter(regions));

      /* std::cout << "num contours: " << regions.size() << std::endl; */
      Contour contour;
      for (const auto& region : regions) {
        orient_contour(region, contour);
        const auto& p = contour.points[0].point();
        const auto& q = contour.points[contour.points.size() - 1].point();
        if (internal::are_equal_points_2(p, q))
          contour.is_closed = true;
        else
          contour.is_closed = false;
        /* std::cout << "is closed: " << contour.is_closed << std::endl; */
        contours.push_back(contour);
      }
    }

    void create_segment_wrappers(
      std::vector<Segment_wrapper>& segs) {
      segs.clear();
      segs.resize(segments.size());

      Indices ns, nt;
      for (std::size_t i = 0; i < segments.size(); ++i) {
        segs[i].index = i;
        segs[i].neighbors.clear();

        const auto& source = segments[i].source();
        const auto& target = segments[i].target();
        find_neighbors(i, source, ns);
        find_neighbors(i, target, nt);

        add_neighbors(ns, segs[i].neighbors);
        add_neighbors(nt, segs[i].neighbors);
      }

      std::sort(segs.begin(), segs.end(),
      [](const Segment_wrapper& a, const Segment_wrapper& b) {
        return a.neighbors.size() < b.neighbors.size();
      });
    }

    void find_neighbors(
      const std::size_t skip,
      const Point_2& query, Indices& neighbors) {
      neighbors.clear();

      for (std::size_t i = 0; i < segments.size(); ++i) {
        if (i == skip) continue;

        const auto& source = segments[i].source();
        const auto& target = segments[i].target();

        if (
          internal::are_equal_points_2(query, source) ||
          internal::are_equal_points_2(query, target) ) {
          neighbors.push_back(i);
        }
      }
    }

    void add_neighbors(
      const Indices& ns, Indices& neighbors) {

      if (ns.size() == 0) return;
      if (ns.size() == 1) {
        neighbors.push_back(ns[0]); return;
      }
    }

    void orient_contour(
      const Indices& region, Contour& contour) {

      contour.points.clear(); My_point mp;
      const std::size_t rs = region.size() - 1;
      for (std::size_t i = 0; i < rs; ++i) {
        const std::size_t ip = i + 1;

        const auto& curr = segments[region[i]];
        const auto& next = segments[region[ip]];

        if (
          internal::are_equal_points_2(curr.source(), next.source()) ||
          internal::are_equal_points_2(curr.source(), next.target()) ) {

          mp.point_ = curr.target();
          mp.labels = curr.lt;
          contour.points.push_back(mp); continue;
        }

        if (
          internal::are_equal_points_2(curr.target(), next.source()) ||
          internal::are_equal_points_2(curr.target(), next.target()) ) {

          mp.point_ = curr.source();
          mp.labels = curr.ls;
          contour.points.push_back(mp); continue;
        }
      }

      const auto& curr = segments[region[rs - 1]];
      const auto& next = segments[region[rs]];

      if (internal::are_equal_points_2(curr.source(), next.source()) ) {

        mp.point_ = curr.source(); mp.labels = curr.ls;
        contour.points.push_back(mp);
        mp.point_ = next.target(); mp.labels = next.lt;
        contour.points.push_back(mp);
        return;
      }

      if (internal::are_equal_points_2(curr.source(), next.target()) ) {

        mp.point_ = curr.source(); mp.labels = curr.ls;
        contour.points.push_back(mp);
        mp.point_ = next.source(); mp.labels = next.ls;
        contour.points.push_back(mp);
        return;
      }

      if (internal::are_equal_points_2(curr.target(), next.source()) ) {

        mp.point_ = curr.target(); mp.labels = curr.lt;
        contour.points.push_back(mp);
        mp.point_ = next.target(); mp.labels = next.lt;
        contour.points.push_back(mp);
        return;
      }

      if (internal::are_equal_points_2(curr.target(), next.target()) ) {

        mp.point_ = curr.target(); mp.labels = curr.lt;
        contour.points.push_back(mp);
        mp.point_ = next.source(); mp.labels = next.ls;
        contour.points.push_back(mp);
        return;
      }
    }

    bool intersect_labels(
      const std::map<std::size_t, Plane_3>& plane_map,
      Line_2& line_2) {

      CGAL_assertion(label_pairs.size() == 1);
      const auto& label_pair = label_pairs[0];

      const auto& plane1 = plane_map.at(label_pair.first);
      const auto& plane2 = plane_map.at(label_pair.second);

      typename CGAL::cpp11::result_of<
      Intersect_3(Plane_3, Plane_3)>::type result
        = CGAL::intersection(plane1, plane2);

      Line_3 line_3; bool found = false;
      if (result) {
        if (const Line_3* l = boost::get<Line_3>(&*result)) {
          found = true; line_3 = *l;
        }
      }
      if (!found) return false;

      const auto p1 = line_3.point(0);
      const auto p2 = line_3.point(1);
      const auto q1 = Point_2(p1.x(), p1.y());
      const auto q2 = Point_2(p2.x(), p2.y());

      line_2 = Line_2(q1, q2);
      return true;
    }

    void set_linear_points(
      const Line_2& line,
      Contour& contour) {

      auto& items = contour.points;
      for (auto& item : items) {

        const auto& p = item.point();
        const auto  q = line.projection(p);

        const FT distance = internal::distance(p, q);
        if (distance < noise_level)
          item.belongs_to_line = true;
      }

      std::size_t start = std::size_t(-1);
      std::size_t end   = std::size_t(-1);

      for (std::size_t i = 1; i < items.size() - 1; ++i) {
        auto& item = items[i];
        if (item.belongs_to_line && start == std::size_t(-1))
          start = i;
        if (item.belongs_to_line)
          end = i;
      }

      if (start == std::size_t(-1) || end == std::size_t(-1))
        return;

      if (start == end) {
        items[start].belongs_to_line = false; return;
      }

      if (start < end && !contour.is_closed) {
        for (std::size_t i = start; i <= end; ++i)
          items[i].belongs_to_line = true;
      }
    }

    void simplify_contour(
      const Line_2& line,
      Contour& contour) {

      Indices free, linear;
      const auto& items = contour.points;
      auto& simplified = contour.simplified;
      simplified.clear();
      for (std::size_t i = 0; i < items.size(); ++i) {

        // Other points.
        if (!items[i].belongs_to_line) {

          // Handle linear.
          if (linear.size() != 0)
            simplify_linear(line, items, linear, simplified);
          free.push_back(i);
          // std::cout << "free: " << i << std::endl;
        } else {

          // Handle free.
          if (free.size() != 0)
            simplify_free(items, free, simplified);
          linear.push_back(i);
          // std::cout << "linear: " << i << std::endl;
        }
      }

      if (free.size() != 0)
        simplify_free(items, free, simplified);
      if (linear.size() != 0)
        simplify_linear(line, items, linear, simplified);

      // Fix unqiue linear points.
      for (std::size_t i = 0; i < simplified.size(); ++i) {
        if (simplified[i].type == Point_type::UNIQUE_LINEAR) {
          add_linear_point(i, line, simplified);
          continue;
        }
      }

      // Fix unique free points.
      for (std::size_t i = 0; i < simplified.size(); ++i) {
        if (simplified[i].type == Point_type::UNIQUE_FREE) {
          if (i > 1) {
            if (simplified[i-1].type == Point_type::LINEAR) {
              add_linear_point(i, line, simplified);
              continue;
            }
          }
          if (i < simplified.size() - 1) {
            if (simplified[i+1].type == Point_type::LINEAR) {
              add_linear_point(i, line, simplified);
              continue;
            }
          }
        }
      }
    }

    void add_linear_point(
      const std::size_t idx,
      const Line_2& line,
      std::vector<My_point>& simplified) {

      auto p = simplified[idx].point();
      simplified[idx].point_ = line.projection(p);
      simplified[idx].type = Point_type::LINEAR;
    }

    void simplify_free(
      const std::vector<My_point>& items,
      Indices& free,
      std::vector<My_point>& simplified) {

      if (free.size() >= 2)
        simplify_polyline(items, free, simplified);
      else {
        if (free.size() == 1) {
          My_point mp;
          mp.point_ = items[free[0]].point();
          mp.type = Point_type::UNIQUE_FREE;
          simplified.push_back(mp);
        }
      }
      free.clear();
    }

    void simplify_linear(
      const Line_2& line,
      const std::vector<My_point>& items,
      Indices& linear,
      std::vector<My_point>& simplified) {

      if (linear.size() >= 2)
        simplify_polyline(line, items, linear, simplified);
      else {
        if (linear.size() == 1) {
          My_point mp;
          mp.point_ = items[linear[0]].point();
          mp.type = Point_type::UNIQUE_LINEAR;
          simplified.push_back(mp);
        }
      }
      linear.clear();
    }

    void simplify_polyline(
      const std::vector<My_point>& items,
      const Indices& polyline,
      std::vector<My_point>& simplified) {

      std::vector<Point_2> points;
      points.reserve(polyline.size());

      for (const std::size_t idx : polyline)
        points.push_back(items[idx].point());

      using Cost = CGAL::Polyline_simplification_2::Squared_distance_cost;
      using Stop = CGAL::Polyline_simplification_2::Stop_above_cost_threshold;

      const double threshold = noise_level / FT(10);

      Cost cost;
      Stop stop(threshold);
      std::vector<Point_2> result;
      CGAL::Polyline_simplification_2::simplify(
        points.begin(), points.end(), cost, stop,
        std::back_inserter(result));

      for (const auto& p : result) {
        My_point mp;
        mp.point_ = p;
        mp.type = Point_type::FREE;
        simplified.push_back(mp);
      }
    }

    void simplify_polyline(
      const Line_2& line,
      const std::vector<My_point>& items,
      const Indices& polyline,
      std::vector<My_point>& simplified) {

      const std::size_t nump = polyline.size();
      const auto& p = items[polyline[0]].point();
      const auto& q = items[polyline[nump - 1]].point();

      const Point_2 s = line.projection(p);
      const Point_2 t = line.projection(q);

      /*
      for (const std::size_t idx : polyline) {
        const auto& pt = items[idx].point();
        My_point mp;
        mp.point_ = pt;
        simplified.push_back(mp);
      } */

      My_point mp;
      mp.point_ = s;
      mp.type = Point_type::LINEAR;
      simplified.push_back(mp);
      mp.point_ = t;
      mp.type = Point_type::LINEAR;
      simplified.push_back(mp);
    }

    void regularize_contour(
      const std::vector<Segment_2>& boundaries,
      const std::vector<Point_pair>& boundary_queries,
      std::shared_ptr<KNQ_pair>& knq_ptr,
      Contour& contour) {

      Regular_segment reg;
      std::vector<Regular_segment> regs;

      const auto& items = contour.simplified;
      regs.reserve(items.size() - 1);

      for (std::size_t i = 0; i < items.size() - 1; ++i) {
        const std::size_t ip = i + 1;

        const auto& p = items[i];
        const auto& q = items[ip];

        add_source_point(p, reg);
        add_target_point(q, reg);
        if (p.type == Point_type::FREE || q.type == Point_type::FREE)
          reg.skip = false;
        else reg.skip = true;
        regs.push_back(reg);
      }

      if (contour.is_closed) {
        regularize_closed_contour(
          boundaries, boundary_queries, knq_ptr, regs);
      }
      else {
        regularize_polyline(
          boundaries, boundary_queries, knq_ptr, regs);
      }

      auto& regularized = contour.regularized;
      regularized.clear();

      My_point mp;
      for (const auto& segment : regs) {
        mp.point_ = segment.source();
        mp.type   = segment.source_type;
        mp.bd_idx = segment.source_bd_idx;

        regularized.push_back(mp);
      }

      const auto& last = regs[regs.size() - 1];
      mp.point_ = last.target();
      mp.type   = last.target_type;
      mp.bd_idx = last.target_bd_idx;

      regularized.push_back(mp);

      /*
      std::cout << "regularized: " << std::endl;
      for (const auto& item : regularized)
        std::cout << int(item.type) << " ";
      std::cout << std::endl;
      for (const auto& item : regularized)
        std::cout << item.bd_idx << " ";
      std::cout << std::endl; */
    }

    void add_source_point(
      const My_point& query, Regular_segment& segment) {

      segment.source_ = query.point();
      segment.source_type = query.type;
      segment.saved_source_type = query.saved_type;
      if (query.type == Point_type::BOUNDARY)
        segment.source_bd_idx = query.bd_idx;
      else
        segment.source_bd_idx = std::size_t(-1);
    }

    void add_target_point(
      const My_point& query, Regular_segment& segment) {

      segment.target_ = query.point();
      segment.target_type = query.type;
      segment.saved_target_type = query.saved_type;
      if (query.type == Point_type::BOUNDARY)
        segment.target_bd_idx = query.bd_idx;
      else
        segment.target_bd_idx = std::size_t(-1);
    }

    void regularize_closed_contour(
      const std::vector<Segment_2>& boundaries,
      const std::vector<Point_pair>& boundary_queries,
      std::shared_ptr<KNQ_pair>& knq_ptr,
      std::vector<Regular_segment>& regs) {

      std::vector< std::vector<Segment_2> > input(1);
      for (const auto& reg : regs)
        input[0].push_back(Segment_2(reg.source(), reg.target()));

      Polygon_regularizer regularizer(
        min_length_2, angle_bound_2, ordinate_bound_2);
      std::vector< std::vector<Seg_pair> > contours;
      create_internal_contours(input, contours);

      std::vector<FT_pair> bounds;
      std::vector<Size_pair> skip;
      std::vector<Segment_2> longest;
      std::vector<Indices> groups;
      regularizer.make_default_groups(
        contours, std::size_t(-1), groups);

      get_multiple_directions(
        true,
        regs, boundaries, boundary_queries, regularizer, knq_ptr,
        contours, bounds, skip, longest, groups);
      regularizer.set_data(bounds, skip, longest, groups);

      regularizer.unify_along_contours(contours);
      regularizer.correct_directions(contours);

      /*
      std::cout << "num longest (close): " <<
        regularizer.get_longest().size() << std::endl;
      for (const std::size_t idx : regularizer.get_groups()[0])
        std::cout << idx << " ";
      std::cout << std::endl; */

      regularizer.regularize_contours(input);

      regs.clear();
      regs.resize(input[0].size());
      for (std::size_t i = 0; i < input[0].size(); ++i) {
        regs[i].source_ = input[0][i].source();
        regs[i].target_ = input[0][i].target();
        regs[i].source_type = Point_type::LINEAR;
        regs[i].target_type = Point_type::LINEAR;
      }
    }

    void create_internal_contours(
      const std::vector< std::vector<Segment_2> >& input,
      std::vector< std::vector<Seg_pair> >& output) {

      output.clear();
      std::vector<Seg_pair> segments;
      for (const auto& contour : input) {
        segments.clear();
        for (const auto& segment : contour)
          segments.push_back(std::make_pair(segment, false));
        output.push_back(segments);
      }
    }

    void get_multiple_directions(
      const bool is_closed,
      const std::vector<Regular_segment>& regs,
      const std::vector<Segment_2>& boundaries,
      const std::vector<Point_pair>& boundary_queries,
      Polygon_regularizer& regularizer,
      std::shared_ptr<KNQ_pair>& knq_ptr,
      std::vector< std::vector<Seg_pair> >& contours,
      std::vector<FT_pair>& bounds,
      std::vector<Size_pair>& skip,
      std::vector<Segment_2>& longest,
      std::vector<Indices>& groups) {

      Segment_2 linear_segment;
      const bool line_found = find_linear_segment(regs, linear_segment);

      if (line_found) {
        longest.push_back(linear_segment);
        bounds.push_back(std::make_pair(FT(45), FT(45)));
        skip.push_back(std::make_pair(std::size_t(-1), std::size_t(-1)));
      }

      std::vector<Segment_2> closest;
      find_closest_segments(
        regs, boundaries, boundary_queries, knq_ptr, closest);

      for (const auto& segment : closest) {
        longest.push_back(segment);
        bounds.push_back(std::make_pair(FT(45), FT(45)));
        skip.push_back(std::make_pair(std::size_t(-1), std::size_t(-1)));
      }

      if (is_closed) {
        create_groups_closed(
          line_found, longest, regs, regularizer, contours, groups);
      } else {
        create_groups_open(
          line_found, longest, regs, regularizer, contours, groups);
      }
    }

    bool find_linear_segment(
      const std::vector<Regular_segment>& regs,
      Segment_2& linear_segment) {

      for (const auto& reg : regs) {
        if (
          reg.saved_source_type == Point_type::LINEAR &&
          reg.saved_target_type == Point_type::LINEAR) {

          const auto& s = reg.source();
          const auto& t = reg.target();
          const auto p = direction.projection(s);
          const auto q = direction.projection(t);

          linear_segment = Segment_2(p, q);
          return true;
        }
      }
      return false;
    }

    void find_closest_segments(
      const std::vector<Regular_segment>& regs,
      const std::vector<Segment_2>& boundaries,
      const std::vector<Point_pair>& boundary_queries,
      std::shared_ptr<KNQ_pair>& knq_ptr,
      std::vector<Segment_2>& closest) {

      Indices neighbors;
      std::set<std::size_t> unique;
      for (const auto& reg : regs) {
        if (reg.source_type == Point_type::BOUNDARY) {
          unique.insert(reg.source_bd_idx); continue;
        }

        const auto& query = reg.source();
        (*knq_ptr)(query, neighbors);
        const std::size_t bd_idx =
          boundary_queries[neighbors[0]].second;
        unique.insert(bd_idx);
      }

      const auto& last = regs[regs.size() - 1];
      if (last.target_type == Point_type::BOUNDARY)
        unique.insert(last.target_bd_idx);

      closest.clear();
      closest.reserve(unique.size());

      for (const std::size_t idx : unique) {
        const auto& segment = boundaries[idx];
        closest.push_back(segment);
      }

      std::sort(closest.begin(), closest.end(),
      [](const Segment_2& a, const Segment_2& b) {
        return a.squared_length() > b.squared_length();
      });
    }

    void create_groups_closed(
      const bool line_found,
      const std::vector<Segment_2>& longest,
      const std::vector<Regular_segment>& regs,
      Polygon_regularizer& regularizer,
      std::vector< std::vector<Seg_pair> >& contours,
      std::vector<Indices>& groups) {

      for (std::size_t i = 0; i < contours[0].size(); ++i) {
        auto& seg_pair = contours[0][i];
        const auto& reg = regs[i];

        if (
          reg.saved_source_type == Point_type::LINEAR &&
          reg.saved_target_type == Point_type::LINEAR) {

          seg_pair.second = true;
          groups[0][i] = std::size_t(-1); continue;
        }

        std::size_t start = 0, end = longest.size();
        if (line_found) end = 1;

        for (std::size_t j = start; j < end; ++j) {
          const FT angle = regularizer.angle_degree_2(
            longest[j], seg_pair.first);
          const FT angle_2 = regularizer.get_angle_2(angle);

          if (
            (CGAL::abs(angle_2) <= regularizer.get_bound_min()) ||
            (CGAL::abs(angle_2) >= regularizer.get_bound_max()) )  {

            seg_pair.second = true;
            groups[0][i] = j; break;
          }
        }
      }
    }

    void create_groups_open(
      const bool line_found,
      const std::vector<Segment_2>& longest,
      const std::vector<Regular_segment>& regs,
      Polygon_regularizer& regularizer,
      std::vector< std::vector<Seg_pair> >& contours,
      std::vector<Indices>& groups) {

      for (std::size_t i = 0; i < contours[0].size(); ++i) {
        auto& seg_pair = contours[0][i];
        const auto& reg = regs[i];

        if (
          reg.saved_source_type == Point_type::LINEAR &&
          reg.saved_target_type == Point_type::LINEAR) {

          seg_pair.second = true;
          groups[0][i] = std::size_t(-1); continue;
        }

        std::size_t start = 0, end = longest.size();
        /* if (line_found && end > 1) start = 1; */

        for (std::size_t j = start; j < end; ++j) {
          const FT angle = regularizer.angle_degree_2(
            longest[j], seg_pair.first);
          const FT angle_2 = regularizer.get_angle_2(angle);

          if (
            (CGAL::abs(angle_2) <= regularizer.get_bound_min()) ||
            (CGAL::abs(angle_2) >= regularizer.get_bound_max()) )  {

            seg_pair.second = true;
            groups[0][i] = j; break;
          }
        }
      }
    }

    void regularize_polyline(
      const std::vector<Segment_2>& boundaries,
      const std::vector<Point_pair>& boundary_queries,
      std::shared_ptr<KNQ_pair>& knq_ptr,
      std::vector<Regular_segment>& regs) {

      std::vector< std::vector<Segment_2> > input(1);
      for (const auto& reg : regs)
        input[0].push_back(Segment_2(reg.source(), reg.target()));

      Polygon_regularizer regularizer(
        min_length_2, angle_bound_2, ordinate_bound_2);
      std::vector< std::vector<Seg_pair> > contours;
      create_internal_contours(input, contours);

      std::vector<FT_pair> bounds;
      std::vector<Size_pair> skip;
      std::vector<Segment_2> longest;
      std::vector<Indices> groups;
      regularizer.make_default_groups(
        contours, std::size_t(-1), groups);

      get_multiple_directions(
        false,
        regs, boundaries, boundary_queries, regularizer, knq_ptr,
        contours, bounds, skip, longest, groups);

      regularizer.unify_along_contours(contours[0], groups[0]);
      regularizer.correct_directions(contours[0], groups[0]);
      regularizer.set_data(bounds, skip, longest, groups);

      /*
      std::cout << "num longest (open): " <<
        regularizer.get_longest().size() << std::endl;
      for (const std::size_t idx : regularizer.get_groups()[0])
        std::cout << idx << " ";
      std::cout << std::endl; */

      regularizer.regularize_polyline(input[0]);

      Regular_segment new_reg;
      std::vector<Regular_segment> new_regs;

      for (std::size_t i = 0; i < input[0].size(); ++i) {
        new_reg.source_ = input[0][i].source();
        new_reg.target_ = input[0][i].target();
        new_reg.source_type = Point_type::LINEAR;
        new_reg.target_type = Point_type::LINEAR;
        new_regs.push_back(new_reg);

        /*
        new_reg.source_ = input[0][i].target();
        new_reg.target_ = input[0][i].target();
        new_regs.push_back(new_reg); */
      }

      new_regs[0].source_type   = regs[0].source_type;
      new_regs[0].source_bd_idx = regs[0].source_bd_idx;

      new_regs[new_regs.size() - 1].target_type   =
        regs[regs.size() - 1].target_type;
      new_regs[new_regs.size() - 1].target_bd_idx =
        regs[regs.size() - 1].target_bd_idx;

      regs.clear(); regs = new_regs;
    }

    void save_original_grid(
      const std::vector<Pixel>& pxs,
      const Color color0,
      const Color color1) {

      std::vector<Point_3> points0, points1;
      points0.reserve(pxs.size());
      points1.reserve(pxs.size());

      for (const auto& px : pxs) {
        const auto& point = px.point;

        if (px.binary == 0)
          points0.push_back(Point_3(point.x(), point.y(), FT(0)));
        else
          points1.push_back(Point_3(point.x(), point.y(), FT(0)));
      }
      Saver saver;
      saver.export_points(points0, color0,
      "/Users/monet/Documents/gf/lod/logs/buildings/tmp/ms_original0");
      saver.export_points(points1, color1,
      "/Users/monet/Documents/gf/lod/logs/buildings/tmp/ms_original1");
    }

    void save_dual_grid(
      const std::vector<Pixel>& pxs,
      const Color color) {

      std::vector<Point_3> points;
      points.reserve(pxs.size());
      for (const auto& px : pxs) {
        const auto& point = px.point;
        points.push_back(Point_3(point.x(), point.y(), FT(0)));
      }
      Saver saver;
      saver.export_points(points, color,
      "/Users/monet/Documents/gf/lod/logs/buildings/tmp/ms_dual");
    }
  };

  using Seed_map             = internal::Seed_property_map;
  using Linear_image_region  = internal::Linear_image_region<Traits, Pixel>;
  using Planar_image_region  = internal::Planar_image_region<Traits, Pixel>;

  Partition_builder_from_image_2(
    const std::vector<Segment_2>& boundary,
    const Triangulation& lod0,
    ImagePointer& image_ptr,
    Partition_2& partition_2,
    const FT noise_level,
    const FT min_length_2,
    const FT angle_bound_2,
    const FT ordinate_bound_2) :
  m_boundary(boundary),
  m_lod0(lod0),
  m_image_ptr(image_ptr),
  m_partition_2(partition_2),
  m_noise_level(noise_level / FT(2)),
  m_min_length_2(min_length_2),
  m_angle_bound_2(angle_bound_2),
  m_ordinate_bound_2(ordinate_bound_2),
  m_pi(static_cast<FT>(CGAL_PI)) {

    m_partition_2.clear();
    create_image();
    create_boundary_knq();
  }

  void build() {

    std::size_t iter = 0;
    do {
      clean_image(); ++iter;
    } while (iter != 2);
    create_label_pairs();
    create_ridges();

    /* auto& ridge = m_ridges[8]; */

    for (auto& ridge : m_ridges)
      ridge.create_contours();
    mark();
    save_original_polylines("original");

    for (auto& ridge : m_ridges)
      ridge.simplify_contours(m_image_ptr->get_plane_map());
    relocate();
    save_simplified_polylines("simplified");

    for (auto& ridge : m_ridges)
      ridge.regularize_contours(
        m_boundary, m_boundary_queries, m_knq_ptr);
    save_regularized_polylines("regularized");
  }

  void create_triangulation() {

    triangulate(
      m_boundary, m_ridges, m_base);
    transform(m_base, m_partition_2);
    save_partition_2(
      "/Users/monet/Documents/gf/lod/logs/buildings/tmp/partition_step_1", false);
  }

  void compute_visibility() {

    apply_visibility(
      m_lod0, m_base);
    transform(m_base, m_partition_2);
    save_partition_2(
      "/Users/monet/Documents/gf/lod/logs/buildings/tmp/partition_step_2", false);
  }

  void label_faces() {

    apply_naive_labeling(
      m_image, m_base);
    correct_labeling(
      m_base);
    transform(m_base, m_partition_2);
    save_partition_2(
      "/Users/monet/Documents/gf/lod/logs/buildings/tmp/partition_step_3", true);
  }

  void save_original_polylines(
    const std::string name) {

    /* const auto& ridge = m_ridges[8]; */

    std::vector<Segment_2> segments;
    for (const auto& ridge : m_ridges) {
      for (const auto& contour : ridge.contours) {
        const auto& items = contour.points;

        for (std::size_t i = 0; i < items.size() - 1; ++i) {
          const std::size_t ip = i + 1;
          const auto& p = items[i].point();
          const auto& q = items[ip].point();
          segments.push_back(Segment_2(p, q));
        }
      }
    }
    Saver saver;
    saver.save_polylines(
      segments, "/Users/monet/Documents/gf/lod/logs/buildings/tmp/contours-" + name);
  }

  void save_simplified_polylines(
    const std::string name) {

    /* const auto& ridge = m_ridges[8]; */

    std::vector<Segment_2> segments;
    for (const auto& ridge : m_ridges) {
      for (const auto& contour : ridge.contours) {
        const auto& items = contour.simplified;

        for (std::size_t i = 0; i < items.size() - 1; ++i) {
          const std::size_t ip = i + 1;
          const auto& p = items[i].point();
          const auto& q = items[ip].point();
          segments.push_back(Segment_2(p, q));
        }
      }
    }
    Saver saver;
    saver.save_polylines(
      segments, "/Users/monet/Documents/gf/lod/logs/buildings/tmp/contours-" + name);
  }

  void save_regularized_polylines(
    const std::string name) {

    /* const auto& ridge = m_ridges[8]; */

    std::vector<Segment_2> segments;
    for (const auto& ridge : m_ridges) {
      for (const auto& contour : ridge.contours) {
        const auto& items = contour.regularized;

        for (std::size_t i = 0; i < items.size() - 1; ++i) {
          const std::size_t ip = i + 1;
          const auto& p = items[i].point();
          const auto& q = items[ip].point();
          segments.push_back(Segment_2(p, q));
        }
      }
    }
    Saver saver;
    saver.save_polylines(
      segments, "/Users/monet/Documents/gf/lod/logs/buildings/tmp/contours-" + name);
  }

  void get_roof_planes(
    std::vector<Plane_3>& roof_planes) {

    const auto& plane_map = m_image_ptr->get_plane_map();
    roof_planes.clear();
    roof_planes.reserve(plane_map.size());

    for (const auto& pair : plane_map) {
      const auto& plane = pair.second;
      roof_planes.push_back(plane);
    }
  }

private:
  const std::vector<Segment_2>& m_boundary;
  const Triangulation& m_lod0;
  ImagePointer& m_image_ptr;
  Partition_2& m_partition_2;
  const FT m_noise_level;
  const FT m_min_length_2;
  const FT m_angle_bound_2;
  const FT m_ordinate_bound_2;
  const FT m_pi;

  Image m_image;
  std::vector<Image> m_ridges;
  std::vector<Point_pair> m_boundary_queries;
  std::shared_ptr<KNQ_pair> m_knq_ptr;
  Triangulation m_base;

  void create_boundary_knq() {

    create_boundary_queries();
    Point_pair_map pmap;
    m_knq_ptr =
      std::make_shared<KNQ_pair>(m_boundary_queries, FT(1), pmap);
  }

  void create_boundary_queries() {

    m_boundary_queries.clear();
    std::vector<Point_2> samples;
    const std::size_t num_samples_per_segment = 20;

    for (std::size_t i = 0; i < m_boundary.size(); ++i) {
      const auto& segment = m_boundary[i];

      const auto& s = segment.source();
      const auto& t = segment.target();

      samples.clear();
      PS_generator generator(s, t, num_samples_per_segment);
      std::copy_n(generator, num_samples_per_segment - 1,
      std::back_inserter(samples));

      for (const auto& p : samples)
        m_boundary_queries.push_back(std::make_pair(p, i));
    }
  }

  void create_image() {

    m_image.clear();
    const auto& original = m_image_ptr->get_image();
    const std::size_t num_labels = m_image_ptr->get_num_labels();
    m_image.num_labels = num_labels;
    auto& pixels = m_image.pixels;
    auto& seeds = m_image.seeds;

    Pixel impixel;
    Idx_map idx_map;

    std::size_t index = 0;
    for (std::size_t i = 0; i < original.rows; ++i) {
      for (std::size_t j = 0; j < original.cols; ++j) {

        const auto& cell = original.grid[i][j];
        const std::size_t label = m_image_ptr->get_label(
          cell.zr, cell.zg, cell.zb);

        if (label == num_labels) {
          impixel.label = std::size_t(-1);
          impixel.is_outer = true;
        } else {
          impixel.label = label;
          impixel.is_outer = false;
        }

        impixel.i = i; impixel.j = j;
        impixel.index = index;
        impixel.point = m_image_ptr->get_point(i, j);
        pixels.push_back(impixel);

        if (impixel.label == std::size_t(-1))
          seeds.push_back(impixel.label);
        else
          seeds.push_back(impixel.index);

        idx_map[std::make_pair(impixel.i, impixel.j)] = impixel.index;
        ++index;
      }
    }

    for (auto& pixel : pixels)
      transform_pixel_neighbors(
        original.rows, original.cols, idx_map, pixel);
  }

  void transform_pixel_neighbors(
    const std::size_t rows, const std::size_t cols,
    const Idx_map& idx_map, Pixel& pixel) {

    get_neighbors_03(
      rows, cols, pixel.i, pixel.j, idx_map, pixel.neighbors_03);
    get_neighbors_47(
      rows, cols, pixel.i, pixel.j, idx_map, pixel.neighbors_47);
  }

  void get_neighbors_03(
    const std::size_t rows, const std::size_t cols,
    const std::size_t i, const std::size_t j,
    const Idx_map& idx_map, Indices& neighbors) {

    neighbors.clear();
    std::size_t ii, jj;

    if (i != 0) {
      ii = i - 1; jj = j;
      const std::size_t idx = idx_map.at(std::make_pair(ii, jj));
      neighbors.push_back(idx);
    }

    if (j != cols - 1) {
      ii = i; jj = j + 1;
      const std::size_t idx = idx_map.at(std::make_pair(ii, jj));
      neighbors.push_back(idx);
    }

    if (i != rows - 1) {
      ii = i + 1; jj = j;
      const std::size_t idx = idx_map.at(std::make_pair(ii, jj));
      neighbors.push_back(idx);
    }

    if (j != 0) {
      ii = i; jj = j - 1;
      const std::size_t idx = idx_map.at(std::make_pair(ii, jj));
      neighbors.push_back(idx);
    }
  }

  void get_neighbors_47(
    const std::size_t rows, const std::size_t cols,
    const std::size_t i, const std::size_t j,
    const Idx_map& idx_map, Indices& neighbors) {

    neighbors.clear();
    std::size_t ii, jj;

    if (i != 0 && j != 0) {
      ii = i - 1;  jj = j - 1;
      const std::size_t idx = idx_map.at(std::make_pair(ii, jj));
      neighbors.push_back(idx);
    }

    if (i != 0 && j != cols - 1) {
      ii = i - 1; jj = j + 1;
      const std::size_t idx = idx_map.at(std::make_pair(ii, jj));
      neighbors.push_back(idx);
    }

    if (i != rows - 1 && j != cols - 1) {
      ii = i + 1; jj = j + 1;
      const std::size_t idx = idx_map.at(std::make_pair(ii, jj));
      neighbors.push_back(idx);
    }

    if (i != rows - 1 && j != 0) {
      ii = i + 1; jj = j - 1;
      const std::size_t idx = idx_map.at(std::make_pair(ii, jj));
      neighbors.push_back(idx);
    }
  }

  void clean_image() {

    using Region_growing = internal::Region_growing<
      Indices, Image, Planar_image_region, Seed_map>;

    auto& pixels = m_image.pixels;
    const auto& seeds = m_image.seeds;

    m_image.use_version_4();
    Seed_map seed_map(seeds);
    Planar_image_region planar_region(pixels);
    Region_growing region_growing(
      seeds, m_image, planar_region, seed_map);

    std::vector<Indices> regions;
    region_growing.detect(std::back_inserter(regions));
    std::cout << "num labeled regions: " << regions.size() << std::endl;

    auto& original = m_image_ptr->get_image();
    for (const auto& region : regions) {
      if (region.size() <= 50) {

        const std::size_t new_label = get_best_label(region);
        if (new_label == std::size_t(-1))
          continue;

        const auto& p = m_image_ptr->get_label_map().at(new_label);
        for (const std::size_t idx : region) {
          const std::size_t i = pixels[idx].i;
          const std::size_t j = pixels[idx].j;

          auto& cell = original.grid[i][j];
          cell.zr = p.x(); cell.zg = p.y(); cell.zb = p.z();
          pixels[idx].label = new_label;
        }
      }
    }

    m_image_ptr->save_image(
      "/Users/monet/Documents/gf/lod/logs/buildings/tmp/image-clean.jpg", original);
  }

  std::size_t get_best_label(
    const Indices& region) {

    const auto& pixels = m_image.pixels;
    const std::size_t num_labels = m_image.num_labels;

    Indices neighbors;
    Indices nums(num_labels, 0);
    for (const std::size_t idx : region) {

      const std::size_t ref_label = pixels[idx].label;
      CGAL_assertion(ref_label != std::size_t(-1));
      if (ref_label == std::size_t(-1))
        continue;

      neighbors.clear();
      m_image(idx, neighbors);

      for (const std::size_t neighbor : neighbors) {
        if (
          pixels[neighbor].label != std::size_t(-1) &&
          pixels[neighbor].label != ref_label) {

          nums[pixels[neighbor].label] += 1;
        }
      }
    }

    std::size_t best_idx = std::size_t(-1);
    std::size_t max_value = 0;
    for (std::size_t i = 0; i < nums.size(); ++i) {
      if (nums[i] > max_value) {
        max_value = nums[i];
        best_idx = i;
      }
    }

    /* CGAL_assertion(best_idx != std::size_t(-1)); */
    return best_idx;
  }

  void create_label_pairs() {

    const auto& pixels = m_image.pixels;
    auto& label_pairs = m_image.label_pairs;

    m_image.use_version_4();
    std::set<Size_pair> unique;
    Indices neighbors;
    for (const auto& pixel : pixels) {
      if (pixel.label == std::size_t(-1))
        continue;

      neighbors.clear();
      m_image(pixel.index, neighbors);

      for (const std::size_t neighbor : neighbors) {
        if (
          pixels[neighbor].label != std::size_t(-1) &&
          pixels[neighbor].label != pixel.label) {

          const std::size_t l1 = pixel.label;
          const std::size_t l2 = pixels[neighbor].label;

          CGAL_assertion(l1 != l2);
          if (l1 <= l2)
            unique.insert(std::make_pair(l1, l2));
          else
            unique.insert(std::make_pair(l2, l1));
        }
      }
    }

    label_pairs.clear();
    label_pairs.reserve(unique.size());
    for (const auto& item : unique)
      label_pairs.push_back(item);
    std::cout << "num label pairs: " << label_pairs.size() << std::endl;
  }

  void create_ridges() {
    m_ridges.clear();
    const auto& label_pairs = m_image.label_pairs;
    for (const auto& label_pair : label_pairs)
      extract_ridges(label_pair);
    std::cout << "num ridges: " << m_ridges.size() << std::endl;
  }

  void extract_ridges(
    const Size_pair& label_pair) {

    Image rimage;
    const bool success = create_ridge_image(
      label_pair, rimage);
    CGAL_assertion(success);
    if (!success) return;
    add_ridges(rimage);
  }

  bool create_ridge_image(
    const Size_pair& label_pair,
    Image& rimage) {

    m_image.use_version_8();
    const auto& pixels = m_image.pixels;
    const std::size_t num_labels = m_image.num_labels;

    Indices neighbors;
    std::set<std::size_t> unique;
    for (const auto& pixel : pixels) {
      if (pixel.label == label_pair.first) {

        neighbors.clear();
        m_image(pixel.index, neighbors);

        bool found = false;
        for (const std::size_t neighbor : neighbors) {
          if (pixels[neighbor].label == label_pair.second) {
            found = true; break;
          }
        }
        if (found) {
          unique.insert(pixel.index);
          for (const std::size_t neighbor : neighbors)
            unique.insert(neighbor);
        }
      }
    }

    if (unique.size() == 0)
      return false;

    std::map<std::size_t, std::size_t> mapping;
    std::size_t index = 0; Pixel rpixel;

    rimage.clear();
    rimage.label_pairs.push_back(label_pair);
    rimage.num_labels = num_labels;

    for (const std::size_t idx : unique) {
      rpixel = pixels[idx];
      mapping[rpixel.index] = index;
      rpixel.index = index;
      ++index;

      rimage.pixels.push_back(rpixel);
      if (rpixel.label == std::size_t(-1))
        rimage.seeds.push_back(rpixel.label);
      else
        rimage.seeds.push_back(rpixel.index);
    }

    for (auto& pixel : rimage.pixels) {
      transform_pixel_neighbors(mapping, pixel.neighbors_03, neighbors);
      pixel.neighbors_03 = neighbors;
      transform_pixel_neighbors(mapping, pixel.neighbors_47, neighbors);
      pixel.neighbors_47 = neighbors;
    }

    /* save_ridge_image(rimage.label_pairs[0], 0, rimage.pixels); */
    return true;
  }

  void transform_pixel_neighbors(
    const std::map<std::size_t, std::size_t>& mapping,
    const Indices& neighbors,
    Indices& result) {

    result.clear();
    for (const std::size_t idx : neighbors) {
      auto neighbor = mapping.find(idx);
      if (neighbor != mapping.end())
        result.push_back(neighbor->second);
      else
        result.push_back(std::size_t(-1));
    }
  }

  void add_ridges(
    Image& rimage) {

    const auto& rpixels = rimage.pixels;
    rimage.use_version_8();
    using Region_growing = internal::Region_growing<
      std::vector<Pixel>, Image, Linear_image_region>;

    Linear_image_region linear_region;
    Region_growing region_growing(
      rpixels, rimage, linear_region);

    std::vector<Indices> regions;
    region_growing.detect(std::back_inserter(regions));

    for (std::size_t i = 0; i < regions.size(); ++i)
      add_ridge(rimage, i, regions);
  }

  void add_ridge(
    const Image& rimage,
    const std::size_t rindex,
    const std::vector<Indices>& regions) {

    const auto& region = regions[rindex];
    const auto& rpixels = rimage.pixels;
    const std::size_t num_labels = rimage.num_labels;
    const auto& label_pairs = rimage.label_pairs;

    Image ridge;
    ridge.noise_level      = m_noise_level;
    ridge.min_length_2     = m_min_length_2;
    ridge.angle_bound_2    = m_angle_bound_2;
    ridge.ordinate_bound_2 = m_ordinate_bound_2;
    ridge.clear();

    ridge.num_labels = num_labels;
    ridge.label_pairs = label_pairs;

    std::map<std::size_t, std::size_t> mapping;
    std::size_t index = 0; Pixel rpixel;

    for (const std::size_t idx : region) {
      rpixel = rpixels[idx];
      mapping[rpixel.index] = index;
      rpixel.index = index;
      ++index;

      ridge.pixels.push_back(rpixel);
      if (rpixel.label == std::size_t(-1))
        ridge.seeds.push_back(rpixel.label);
      else
        ridge.seeds.push_back(rpixel.index);
    }

    Indices neighbors;
    for (auto& pixel : ridge.pixels) {
      transform_pixel_neighbors(mapping, pixel.neighbors_03, neighbors);
      pixel.neighbors_03 = neighbors;
      transform_pixel_neighbors(mapping, pixel.neighbors_47, neighbors);
      pixel.neighbors_47 = neighbors;
    }

    /* save_ridge_image(ridge.label_pairs[0], rindex, ridge.pixels); */

    ridge.is_ridge = true;
    m_ridges.push_back(ridge);
  }

  void mark() {

    for (std::size_t i = 0; i < m_ridges.size(); ++i) {
      for (auto& contour : m_ridges[i].contours) {
        if (contour.is_closed) continue;
        const std::size_t nump = contour.points.size();
        auto& neighbors = contour.neighbors;

        for (std::size_t k = 0; k < nump; ++k) {
          const auto& query = contour.points[k];
          search_for_neighbors(i, query, neighbors);
        }
      }
    }

    std::vector<Size_pair> ns;
    for (std::size_t i = 0; i < m_ridges.size(); ++i) {
      for (auto& contour : m_ridges[i].contours) {
        if (contour.is_closed) continue;
        const auto& neighbors = contour.neighbors;

        auto& items = contour.points;
        const std::size_t nump = items.size();

        for (std::size_t k = 0; k < nump; ++k) {
          const auto& query = items[k];
          find_corner_neighbors(query, neighbors, ns);

          if (k > nump / 2 - 1) {
            for (const auto& n : ns)
              contour.ns1.insert(n);
          } else {
            for (const auto& n : ns)
              contour.ns0.insert(n);
          }
        }
      }
    }
  }

  void search_for_neighbors(
    const std::size_t skip, const My_point& query,
    std::set<Size_pair>& neighbors) {

    for (std::size_t i = 0; i < m_ridges.size(); ++i) {
      if (i == skip) continue;

      for (std::size_t j = 0; j < m_ridges[i].contours.size(); ++j) {
        const auto& contour = m_ridges[i].contours[j];

        if (contour.is_closed) continue;
        const std::size_t nump = contour.points.size();

        for (std::size_t k = 0; k < nump; ++k) {
          const auto& item = contour.points[k];
          if (internal::are_equal_points_2(item.point(), query.point()))
            neighbors.insert(std::make_pair(i, j));
        }
      }
    }
  }

  void find_corner_neighbors(
    const My_point& query,
    const std::set<Size_pair>& neighbors,
    std::vector<Size_pair>& ns) {

    ns.clear();
    for (const auto& neighbor : neighbors) {
      const auto& contour = m_ridges[neighbor.first].contours[neighbor.second];
      const auto& items = contour.points;

      for (const auto& item : items) {
        if (internal::are_equal_points_2(query.point(), item.point())) {
          ns.push_back(neighbor); break;
        }
      }
    }
  }

  void relocate() {
    apply_point_types();
    relocate_corners();
    relocate_boundary_points();
  }

  void apply_point_types() {

    for (auto& ridge : m_ridges) {
      for (auto& contour : ridge.contours) {
        auto& items = contour.simplified;
        for (auto& item : items)
          item.saved_type = item.type;
      }
    }

    for (auto& ridge : m_ridges) {
      for (auto& contour : ridge.contours) {
        if (contour.is_closed) continue;

        auto& items = contour.simplified;
        const std::size_t nump = items.size();
        auto& p = items[0];
        auto& q = items[nump - 1];

        if (contour.ns0.size() == 0)
          p.type = Point_type::BOUNDARY;
        else
          p.type = Point_type::CORNER;

        if (contour.ns1.size() == 0)
          q.type = Point_type::BOUNDARY;
        else
          q.type = Point_type::CORNER;
      }
    }
  }

  void relocate_corners() {
    for (auto& ridge : m_ridges) {
      for (auto& contour : ridge.contours) {
        if (contour.is_closed) continue;

        auto& items = contour.simplified;
        const std::size_t nump = items.size();
        auto& p = items[0];
        auto& q = items[nump - 1];

        if (p.type == Point_type::CORNER)
          relocate_corner(contour.ns0, p);
        if (q.type == Point_type::CORNER)
          relocate_corner(contour.ns1, q);
      }
    }
  }

  void relocate_boundary_points() {
    for (auto& ridge : m_ridges) {
      for (auto& contour : ridge.contours) {
        if (contour.is_closed) continue;

        auto& items = contour.simplified;
        const std::size_t nump = items.size();
        auto& p = items[0];
        auto& q = items[nump - 1];

        if (p.type == Point_type::BOUNDARY)
          intersect_with_boundary(p);
        if (q.type == Point_type::BOUNDARY)
          intersect_with_boundary(q);
      }
    }
  }

  void intersect_with_boundary(
    My_point& query) {

    Indices closest;
    (*m_knq_ptr)(query.point(), closest);
    const std::size_t bd_idx = m_boundary_queries[closest[0]].second;

    const auto& segment = m_boundary[bd_idx];
    const auto& s = segment.source();
    const auto& t = segment.target();
    const Line_2 line2 = Line_2(s, t);

    query.bd_idx = bd_idx;
    query.point_ = line2.projection(query.point());
  }

  bool intersect_2(
    const Line_2& line_1, const Line_2& line_2,
    Point_2& in_point) {

    typename std::result_of<Intersect_2(Line_2, Line_2)>::type result
    = CGAL::intersection(line_1, line_2);
    if (result) {
      if (const Line_2* line = boost::get<Line_2>(&*result))
        return false;
      else {
        const Point_2* point = boost::get<Point_2>(&*result);
        in_point = *point; return true;
      }
    }
    return false;
  }

  void relocate_corner(
    const std::set<Size_pair>& ns,
    My_point& query) {

    Point_2 new_pos;
    compute_centered_corner(query, ns, new_pos);
    apply_new_corner_position(ns, new_pos, query);
  }

  void compute_centered_corner(
    const My_point& query,
    const std::set<Size_pair>& ns,
    Point_2& center) {

    FT x = query.point().x();
    FT y = query.point().y();

    for (const auto& n : ns) {
      const auto& contour = m_ridges[n.first].contours[n.second];
      const auto& items = contour.simplified;
      const std::size_t bd_idx = get_bd_idx(query, items);

      const auto& other = items[bd_idx];
      x += other.point().x();
      y += other.point().y();
    }

    x /= static_cast<FT>(ns.size() + 1);
    y /= static_cast<FT>(ns.size() + 1);

    center = Point_2(x, y);
  }

  std::size_t get_bd_idx(
    const My_point& query,
    const std::vector<My_point>& items) {

    const std::size_t nump = items.size();

    const auto& p = items[0];
    const auto& q = items[nump - 1];

    const FT dist1 = internal::distance(query.point(), p.point());
    const FT dist2 = internal::distance(query.point(), q.point());

    if (dist1 < dist2) return 0;
    return nump - 1;
  }

  void apply_new_corner_position(
    const std::set<Size_pair>& ns,
    const Point_2& new_pos,
    My_point& query) {

    for (const auto& n : ns) {
      auto& contour = m_ridges[n.first].contours[n.second];
      auto& items = contour.simplified;
      const std::size_t bd_idx = get_bd_idx(query, items);

      auto& other = items[bd_idx];
      other.point_ = new_pos;
    }
    query.point_ = new_pos;
  }

  void triangulate(
    const std::vector<Segment_2>& boundary,
    const std::vector<Image>& ridges,
    Triangulation& base) {

    auto& tri = base.delaunay;
    tri.clear();

    for (const auto& segment : boundary) {
      const auto vh1 = tri.insert(segment.source());
      const auto vh2 = tri.insert(segment.target());

      if (vh1 != vh2)
        tri.insert_constraint(vh1, vh2);
    }

    std::vector<Vertex_handle> vhs;
    for (const auto& ridge : ridges) {

      vhs.clear();
      for (const auto& contour : ridge.contours) {

        const auto& items = contour.regularized;
        for (const auto& item : items)
          vhs.push_back(tri.insert(item.point()));

        for (std::size_t i = 0; i < items.size() - 1; ++i) {
          const std::size_t ip = i + 1;
          if (vhs[i] != vhs[ip])
            tri.insert_constraint(vhs[i], vhs[ip]);
        }
      }
    }
  }

  void apply_visibility(
    const Triangulation& lod0,
    Triangulation& base) {

    for (auto fh = base.delaunay.finite_faces_begin();
    fh != base.delaunay.finite_faces_end(); ++fh) {

      const Point_2 center = CGAL::barycenter(
        fh->vertex(0)->point(), FT(1),
        fh->vertex(1)->point(), FT(1),
        fh->vertex(2)->point(), FT(1));

      const auto handle = lod0.delaunay.locate(center);
      if (handle->info().tagged)
        fh->info().interior = true;
      else
        fh->info().interior = false;
    }
  }

  void apply_naive_labeling(
    const Image& image,
    Triangulation& base) {

    const std::size_t num_labels = m_image_ptr->get_num_labels();
    auto& tri = base.delaunay;
    for (auto fh = tri.finite_faces_begin();
    fh != tri.finite_faces_end(); ++fh) {
      if (fh->info().interior) {

        fh->info().probabilities.clear();
        fh->info().probabilities.resize(num_labels, FT(0));
      }
    }

    for (const auto& pixel : image.pixels) {
      if (pixel.label == std::size_t(-1))
        continue;

      const auto& p = pixel.point;
      const auto fh = tri.locate(p);
      if (fh->info().interior)
        fh->info().probabilities[pixel.label] += FT(1);
    }

    for (auto fh = tri.finite_faces_begin();
    fh != tri.finite_faces_end(); ++fh) {
      if (fh->info().interior) {

        FT max_prob = FT(-1);
        std::size_t best_label = std::size_t(-1);
        for (std::size_t i = 0; i < num_labels; ++i) {
          const FT prob = fh->info().probabilities[i];
          if (prob > max_prob) {
            max_prob = prob;
            best_label = i;
          }
        }

        if (max_prob != FT(0))
          fh->info().label = best_label;
        else
          fh->info().label = std::size_t(-1);
      }
    }
  }

  void correct_labeling(
    Triangulation& base) {

    apply_ray_shooting(base);
    fill_holes(base);
    recolor_unique(base);
  }

  void apply_ray_shooting(
    Triangulation& base) {

    const FT radius = FT(1);
    const std::size_t num_samples = 24;

    std::vector<Point_2> samples1, samples2;
    samples1.reserve(num_samples / 2);
    samples2.reserve(num_samples / 2);

    std::vector<Point_2> samples;
    samples.reserve(num_samples);

    auto& tri = base.delaunay;
    for (auto fh = tri.finite_faces_begin();
    fh != tri.finite_faces_end(); ++fh) {
      if (fh->info().interior) {

        const Point_2 p = CGAL::barycenter(
          fh->vertex(0)->point(), FT(1),
          fh->vertex(1)->point(), FT(1),
          fh->vertex(2)->point(), FT(1));

        if (tri.oriented_side(fh, p) == CGAL::ON_NEGATIVE_SIDE)
          continue;

        create_points_on_circle(
          p, radius, FT(0), num_samples, samples1);
        create_points_on_circle(
          p, radius, FT(180), num_samples, samples2);

        samples.clear();
        for (const auto& sample : samples1)
          samples.push_back(sample);
        for (const auto& sample : samples2)
          samples.push_back(sample);

        for (const auto& q : samples) {
          LF_circulator circ = tri.line_walk(p, q, fh);
          const LF_circulator end = circ;
          if (circ.is_empty()) continue;

          std::size_t curr = std::size_t(-1);
          std::size_t count = 0;
          do {
            LF_circulator f1 = circ; ++circ;
            LF_circulator f2 = circ;

            if (count != 0) {
              if (f2->info().label != curr)
                break;
            }

            const bool success = are_neighbors(f1, f2);
            if (!success) break;

            const std::size_t idx = f1->index(f2);
            const auto edge = std::make_pair(f1, idx);

            if (tri.is_constrained(edge)) break;
            if (tri.is_infinite(f2)) break;

            if (f2->info().label != std::size_t(-1)) {

              curr = f2->info().label; ++count;
              fh->info().probabilities[curr] += FT(1);
            }

          } while (circ != end);
        }

        FT max_prob = FT(-1);
        std::size_t best_label = std::size_t(-1);
        const std::size_t num_labels = fh->info().probabilities.size();

        for (std::size_t i = 0; i < num_labels; ++i) {
          const FT prob = fh->info().probabilities[i];
          if (prob > max_prob) {
            max_prob = prob; best_label = i;
          }
        }
        if (max_prob != FT(0))
          fh->info().label = best_label;
      }
    }
  }

  void fill_holes(
    Triangulation& base) {

    auto& tri = base.delaunay;
    for (auto fh = tri.finite_faces_begin();
    fh != tri.finite_faces_end(); ++fh) {

      if (fh->info().interior &&
      fh->info().label == std::size_t(-1)) {
        const std::size_t num_labels = fh->info().probabilities.size();

        std::vector<std::size_t> nums(num_labels, 0);
        for (std::size_t k = 0; k < 3; ++k)
          if (!tri.is_infinite(fh->neighbor(k)) &&
            fh->neighbor(k)->info().label != std::size_t(-1))
          nums[fh->neighbor(k)->info().label] += 1;

        std::size_t max_val = 0;
        std::size_t best_label = std::size_t(-1);

        for (std::size_t i = 0; i < num_labels; ++i) {
          const FT val = nums[i];
          if (val > max_val) {
            max_val = val; best_label = i;
          }
        }

        /* CGAL_assertion(best_label != std::size_t(-1)); */

        if (best_label != std::size_t(-1))
          fh->info().label = best_label;
        else
          fh->info().label = 0;
      }
    }
  }

  void recolor_unique(
    Triangulation& base) {

    auto& tri = base.delaunay;
    for (auto fh = tri.finite_faces_begin();
    fh != tri.finite_faces_end(); ++fh) {
      if (fh->info().interior) {

        std::set<std::size_t> unique;
        for (std::size_t k = 0; k < 3; ++k)
          if (!tri.is_infinite(fh->neighbor(k)))
            unique.insert(fh->neighbor(k)->info().label);
        if (unique.size() == 1)
          fh->info().label = *(unique.begin());

        bool found = false;
        for (const auto& val : unique) {
          if (val == fh->info().label) {
            found = true; break;
          }
        }

        if (!found) {
          fh->info().label = *(unique.begin());
        }
      }
    }
  }

  void create_points_on_circle(
    const Point_2& center,
    const FT radius,
    const FT start,
    const std::size_t num_samples,
    std::vector<Point_2>& samples) {

    samples.clear();
    FT factor = FT(360) / static_cast<FT>(num_samples);
    factor *= m_pi; factor /= FT(180);

    FT init = start;
    init *= m_pi; init /= FT(180);

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

  bool are_neighbors(
    LF_circulator f1, LF_circulator f2) const {

    for (std::size_t i = 0; i < 3; ++i) {
      const std::size_t ip = (i + 1) % 3;

      const auto p1 = f1->vertex(i);
      const auto p2 = f1->vertex(ip);

      for (std::size_t j = 0; j < 3; ++j) {
        const std::size_t jp = (j + 1) % 3;

        const auto q1 = f2->vertex(j);
        const auto q2 = f2->vertex(jp);

        if (
          ( p1 == q1 && p2 == q2) ||
          ( p1 == q2 && p2 == q1) ) {

          return true;
        }
      }
    }
    return false;
  }

  void transform(
    const Triangulation& base,
    Partition_2& partition_2) {

    partition_2.clear();
    if (base.delaunay.number_of_faces() < 1)
      return;

    std::map<Face_handle, int> fmap;
    create_faces(base, partition_2, fmap);
    create_face_neighbors(base, fmap, partition_2);
    create_edges(base, fmap, partition_2);
  }

  void create_faces(
    const Triangulation& base,
    Partition_2& partition_2,
    std::map<Face_handle, int>& fmap) const {

    const auto& tri = base.delaunay;

    partition_2.faces.clear();
    partition_2.faces.reserve(tri.number_of_faces());

    fmap.clear();

    Partition_face_2 pface;
    std::vector<Vertex_handle> vhs(3);

    int idx = 0;
    for (auto fit = tri.finite_faces_begin();
    fit != tri.finite_faces_end(); ++fit) {
      const Face_handle fh = static_cast<Face_handle>(fit);
      pface.base.delaunay.clear();

      for (std::size_t k = 0; k < 3; ++k)
        vhs[k] = pface.base.delaunay.insert(fh->vertex(k)->point());
      for (std::size_t k = 0; k < 3; ++k) {
        const std::size_t kp = (k + 1) % 3;
        if (vhs[k] != vhs[kp])
          pface.base.delaunay.insert_constraint(vhs[k], vhs[kp]);
      }

      if (fh->info().interior) {
        pface.inside = FT(1); pface.outside = FT(0);
        pface.visibility = Visibility_label::INSIDE;
        pface.label = fh->info().label;
      } else {
        pface.inside = FT(0); pface.outside = FT(1);
        pface.visibility = Visibility_label::OUTSIDE;
        pface.label = std::size_t(-1);
      }

      partition_2.faces.push_back(pface);
      fmap[fh] = idx; ++idx;
    }
  }

  void create_face_neighbors(
    const Triangulation& base,
    const std::map<Face_handle, int>& fmap,
    Partition_2& partition_2) const {

    const auto& tri = base.delaunay;

    int idx = 0;
    for (auto fit = tri.finite_faces_begin();
    fit != tri.finite_faces_end(); ++fit) {
      const Face_handle fh = static_cast<Face_handle>(fit);

      auto& edges = partition_2.faces[idx].edges;
      auto& neighbors = partition_2.faces[idx].neighbors;

      edges.clear(); edges.reserve(3);
      neighbors.clear(); neighbors.reserve(3);

      for (std::size_t k = 0; k < 3; ++k) {
        const Face_handle fhn = fh->neighbor(k);
        if (tri.is_infinite(fhn)) neighbors.push_back(-1);
        else {
          CGAL_assertion(fmap.find(fhn) != fmap.end());
          neighbors.push_back(fmap.at(fhn));
        }

        const auto& p1 = fh->vertex((k + 1) % 3)->point();
        const auto& p2 = fh->vertex((k + 2) % 3)->point();
        edges.push_back(Segment_2(p1, p2));
      }
      ++idx;
    }
  }

  void create_edges(
    const Triangulation& base,
    const std::map<Face_handle, int>& fmap,
    Partition_2& partition_2) const {

    const auto& tri = base.delaunay;

    partition_2.edges.clear();
    partition_2.edges.reserve(tri.number_of_faces());
    for (auto eh = tri.finite_edges_begin();
    eh != tri.finite_edges_end(); ++eh) {
      const Face_handle fh = eh->first;
      const std::size_t idx = eh->second;
      const Face_handle fhn = fh->neighbor(idx);

      const auto& p1 = fh->vertex((idx + 1) % 3)->point();
      const auto& p2 = fh->vertex((idx + 2) % 3)->point();

      int f1 = -1, f2 = -1;
      if (!tri.is_infinite(fh)) {
        if (fmap.find(fh) != fmap.end())
          f1 = fmap.at(fh);
      }
      if (!tri.is_infinite(fhn)) {
        if (fmap.find(fhn) != fmap.end())
          f2 = fmap.at(fhn);
      }
      partition_2.edges.push_back(Partition_edge_2(p1, p2, f1, f2));
    }
  }

  void save_ridge_image(
    const Size_pair& label_pair,
    const std::size_t ridge_index,
    const std::vector<Pixel>& pixels) {

    const auto& original = m_image_ptr->get_image();
    auto image = original;
    for (const auto& pixel : pixels) {
      auto& cell = image.grid[pixel.i][pixel.j];
      cell.zr = FT(0); cell.zg = FT(0); cell.zb = FT(0);
    }
    m_image_ptr->save_image(
      "/Users/monet/Documents/gf/lod/logs/buildings/tmp/ridges/ridge-"
      + std::to_string(label_pair.first)  + "-"
      + std::to_string(label_pair.second) + "-"
      + std::to_string(ridge_index) + ".jpg", image);
  }

  void save_partition_2(
    const std::string path,
    const bool with_roof_colors) {

    const FT z = FT(0);
    std::size_t num_vertices = 0;
    internal::Indexer<Point_3> indexer;

    std::vector<Point_3> vertices;
    std::vector<Indices> faces;
    std::vector<CGAL::Color> fcolors;

    Polygon_inserter<Traits> inserter(faces, fcolors);
    auto output_vertices = std::back_inserter(vertices);
    auto output_faces = boost::make_function_output_iterator(inserter);

    for (const auto& face : m_partition_2.faces) {
      if (!with_roof_colors) {
        face.output_for_visibility(
          indexer, num_vertices, output_vertices, output_faces, z);
      } else {
        if (face.visibility == Visibility_label::INSIDE)
          face.output_with_label_color(
            indexer, num_vertices, output_vertices, output_faces, z);
      }
    }

    Saver saver;
    saver.export_polygon_soup(vertices, faces, fcolors, path);
  }
};

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_INTERNAL_PARTITION_BUILDER_FROM_IMAGE_2_H
