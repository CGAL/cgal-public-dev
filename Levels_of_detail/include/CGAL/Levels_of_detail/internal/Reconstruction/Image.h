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

#ifndef CGAL_LEVELS_OF_DETAIL_INTERNAL_IMAGE_H
#define CGAL_LEVELS_OF_DETAIL_INTERNAL_IMAGE_H

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

// Internal includes.
#include <CGAL/Levels_of_detail/internal/utils.h>
#include <CGAL/Levels_of_detail/internal/struct.h>

// Other includes.
#include <CGAL/Levels_of_detail/internal/Shape_detection/Linear_image_region.h>
#include <CGAL/Levels_of_detail/internal/Shape_detection/Oriented_image_region.h>
#include <CGAL/Levels_of_detail/internal/Spatial_search/Oriented_neighbor_query.h>

// Testing.
#include "../../../../../test/Levels_of_detail/include/Saver.h"

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

template<
typename GeomTraits>
struct Image {

public:
  using Traits = GeomTraits;

  using FT = typename Traits::FT;
  using Point_2 = typename Traits::Point_2;
  using Point_3 = typename Traits::Point_3;
  using Line_2 = typename Traits::Line_2;
  using Line_3 = typename Traits::Line_3;
  using Segment_2 = typename Traits::Segment_2;
  using Plane_3 = typename Traits::Plane_3;
  using Intersect_3 = typename Traits::Intersect_3;

  using Color = CGAL::Color;

  using Size_pair = std::pair<std::size_t, std::size_t>;
  using Indices = std::vector<std::size_t>;
  using Label_pairs = std::vector<Size_pair>;

  using Saver = Saver<Traits>;

  struct Pixel {
    Point_2 point;
    Indices neighbors_03;
    Indices neighbors_47;
    bool is_outer = false;
    std::size_t i = std::size_t(-1), j = std::size_t(-1);
    std::size_t index = std::size_t(-1);
    std::size_t label = std::size_t(-1);
    std::size_t binary = std::size_t(-1);
    bool used = false;
  };

  using Pixels = std::vector<Pixel>;

  struct Pixel_point_map {

    using key_type = const Pixel&;
    using value_type = Point_2;
    using reference = const value_type&;
    using category = boost::lvalue_property_map_tag;

    friend inline reference get(
      const Pixel_point_map& pmap, key_type pixel) {
      return pixel.point;
    }
  };

  Pixel_point_map pixel_point_map;

  enum class Point_type {
    DEFAULT = 0,
    UNIQUE_FREE = 1,
    UNIQUE_LINEAR = 2,
    FREE = 3,
    LINEAR = 4,
    CORNER = 5,
    BOUNDARY = 6,
    DISTANT_BOUNDARY = 7,
    OUTER_CORNER = 8,
    OUTER_BOUNDARY = 9
  };

  struct My_point {
    Point_2 point_;
    std::set<std::size_t> labels;
    bool belongs_to_line = false;
    Point_type int_type = Point_type::DEFAULT;
    Point_type end_type = Point_type::DEFAULT;
    std::size_t bd_idx = std::size_t(-1);
    std::set<Size_pair> neighbors;
    bool is_removed = false;

    const Point_2& point() const {
      return point_;
    }

    void merge(const My_point& other) {
      const auto mid = internal::middle_point_2(point_, other.point_);
      point_ = mid;

      std::set<std::size_t> labs;
      for (const std::size_t lab : labels)
        labs.insert(lab);
      for (const std::size_t lab : other.labels)
        labs.insert(lab);
      labels = labs;

      std::set<Size_pair> neighs;
      for (const auto& neigh : neighbors)
        neighs.insert(neigh);
      for (const auto& neigh : other.neighbors)
        neighs.insert(neigh);
      neighbors = neighs;
    }

    void add_neighbor(
      const std::size_t i, const std::size_t j) {
      neighbors.insert(std::make_pair(i, j));
    }
  };

  using Points = std::vector<My_point>;

  struct Contour {
    Points points;
    bool is_closed = false;
    bool is_degenerated = false;

    bool skip() const {
      return is_closed || is_degenerated;
    }

    void truncate() {

      const std::size_t nump = points.size();
      const std::size_t mid = nump / 2 - 1;

      std::vector<bool> to_be_removed(nump, false);
      for (std::size_t i = 0; i < nump - 1; ++i) {

        auto& curr = points[i];
        auto& next = points[i + 1];

        auto& ns0 = curr.neighbors;
        auto& ns1 = next.neighbors;

        if (i <= mid) {
          if (ns0.size() >= 2 && ns1.size() < 2) {
            for (std::size_t k = 0; k < i; ++k)
              to_be_removed[k] = true;
            continue;
          } /* else ns0.clear(); */
        } else {
          if (ns0.size() >= 2) {
            for (std::size_t k = i + 1; k < nump; ++k)
              to_be_removed[k] = true;
          } /* else ns0.clear(); */
        }
      }

      /*
      if (points[nump - 1].neighbors.size() < 2)
        points[nump - 1].neighbors.clear(); */

      Points truncated;
      for (std::size_t i = 0; i < to_be_removed.size(); ++i)
        if (!to_be_removed[i])
          truncated.push_back(points[i]);
      points.clear();
      points = truncated;

      /*
      for (const auto& p : points)
        std::cout << p.neighbors.size() << " ";
      std::cout << std::endl << std::endl; */
    }
  };

  using Contours = std::vector<Contour>;

  struct My_segment {
    Point_2 source_;
    Point_2 target_;

    const Point_2& source() const {
      return source_;
    }
    const Point_2& target() const {
      return target_;
    }
  };

  using Segments = std::vector<My_segment>;

  struct Segment_wrapper {
    std::size_t index = std::size_t(-1);
    Indices neighbors;
  };

  using Segment_wrappers = std::vector<Segment_wrapper>;

  using Seed_map =
    internal::Seed_property_map;
  using Linear_image_region =
    internal::Linear_image_region<Traits, Pixel>;
  using Oriented_neighbor_query =
    internal::Oriented_neighbor_query<Traits, Segment_wrapper>;
  using Oriented_image_region =
    internal::Oriented_image_region<Traits, Segment_wrapper>;

  // Variables.
  Indices seeds;
  Pixels pixels;
  Label_pairs label_pairs;
  Contours contours;
  Pixels dual;
  Segments segments;
  Line_2 line;
  bool line_found;
  std::size_t num_labels;
  bool is_ridge;
  FT noise_level_2;

  // Functions.
  void clear() {
    seeds.clear();
    pixels.clear();
    label_pairs.clear();
    contours.clear();
    dual.clear();
    segments.clear();
    line = Line_2();
    line_found = false;
    num_labels = 0;
    is_ridge = false;
    noise_level_2 = FT(-1);
  }

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

  void create_line(
    const std::map<std::size_t, Plane_3>& plane_map) {

    CGAL_assertion(label_pairs.size() == 1);
    const std::size_t l1 = label_pairs[0].first;
    const std::size_t l2 = label_pairs[0].second;
    if (l1 == std::size_t(-1) || l2 == std::size_t(-1)) {
      line_found = false; return;
    }
    line_found = intersect_labels(plane_map, l1, l2, line);
  }

  void create_contours() {

    if (!is_ridge) return;

    make_binary_indices();
    make_dual_grid();
    apply_contouring();
    make_contours();
    type_contours();
    set_labels();
  }

  void create_contours_from_label_pair_v1(
    const std::size_t lp_index,
    const Size_pair& label_pair,
    std::vector<Contour>& result) {

    /* if (lp_index != 10) return; */
    for (auto& pixel : pixels) {
      pixel.binary = std::size_t(-1);
      pixel.used = false;
    }

    make_binary_indices(label_pair);
    make_dual_grid(lp_index, label_pair);
    apply_contouring(lp_index);
    make_contours();
    type_contours();
    set_labels(label_pair);
    result = contours;
  }

  void create_contours_from_label_pair_v2(
    const std::size_t lp_index,
    const Size_pair& label_pair,
    std::vector<Contour>& result) {

    for (auto& pixel : pixels) {
      pixel.binary = std::size_t(-1);
      pixel.used = false;
    }

    const std::size_t l1 = label_pair.first;
    const std::size_t l2 = label_pair.second;

    segments.clear();
    My_segment my_segment;
    for (const auto& pixel1 : pixels) {
      if (pixel1.label != l1) continue;

      const auto& neighbors = pixel1.neighbors_03;
      for (std::size_t i = 0; i < neighbors.size(); ++i) {
        const std::size_t neighbor = neighbors[i];

        if (neighbor == std::size_t(-1)) continue;
        const auto& pixel2 = pixels[neighbor];
        if (pixel2.label != l2) continue;

        auto segment = Segment_2(pixel1.point, pixel2.point);
        rotate(FT(90), segment);
        my_segment.source_ = segment.source();
        my_segment.target_ = segment.target();
        segments.push_back(my_segment);
      }
    }

    /*
    std::vector<Segment_2> contour;
    contour.reserve(segments.size());
    for (const auto& segment : segments)
      contour.push_back(Segment_2(segment.source(), segment.target()));

    Saver saver;
    saver.save_polylines(
      contour, "/Users/monet/Documents/gf/lod/logs/buildings/tmp/duals/ms_contour-" +
      std::to_string(lp_index)); */

    make_contours();
    type_contours();
    set_labels(label_pair);
    result = contours;
  }

private:
  bool m_use_version_8 = false;

  void rotate(
    const FT angle_deg,
    Segment_2& segment) {

    Point_2 source = segment.source();
    Point_2 target = segment.target();

    const FT pi = static_cast<FT>(CGAL_PI);
    const Point_2 b = internal::middle_point_2(source, target);
    const FT angle_rad = angle_deg * pi / FT(180);

    internal::rotate_point_2(angle_rad, b, source);
    internal::rotate_point_2(angle_rad, b, target);

    segment = Segment_2(source, target);
  }

  bool intersect_labels(
    const std::map<std::size_t, Plane_3>& plane_map,
    const std::size_t l1, const std::size_t l2,
    Line_2& line_2) {

    const auto& plane1 = plane_map.at(l1);
    const auto& plane2 = plane_map.at(l2);

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

  void make_binary_indices(const Size_pair& label_pair) {

    const std::size_t l1 = label_pair.first;
    const std::size_t l2 = label_pair.second;

    for (auto& pixel : pixels) {
      if (pixel.label == l1)
        pixel.binary = 0;
      if (pixel.label == l2)
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

    /* save_dual_grid(
      dual, Color(0, 125, 0)); */
  }

  void make_dual_grid(
    const std::size_t lp_index,
    const Size_pair& label_pair) {

    /*
    save_original_grid(pixels,
      Color(125, 0, 0), Color(0, 0, 125), lp_index); */

    std::vector<Point_2> points;
    for (auto& pixel1 : pixels) {
      if (pixel1.binary == std::size_t(-1)) continue;

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
      bool found = false;
      CGAL_assertion(neighbors47.size() == 4);
      for (std::size_t i = 0; i < neighbors47.size(); ++i) {
        const std::size_t neighbor = neighbors47[i];
        if (neighbor == std::size_t(-1)) continue;

        const auto& pixel2 = pixels[neighbor];
        if (pixel2.used) continue;
        if (pixel2.binary == std::size_t(-1)) continue;
        if (pixel2.binary != pixel1.binary)
          found = true;
      }
      if (!found) continue;

      for (std::size_t i = 0; i < neighbors47.size(); ++i) {
        const std::size_t neighbor = neighbors47[i];
        if (neighbor == std::size_t(-1)) continue;
        const auto& pixel2 = pixels[neighbor];

        if (pixel2.used) continue;
        if (pixel2.binary == std::size_t(-1)) continue;

        const auto& p = pixel1.point;
        const auto& q = pixel2.point;
        const auto  m = internal::middle_point_2(p, q);
        points.push_back(m);
      }
    }

    /* std::cout << "points: " << points.size() << std::endl; */
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

    std::size_t count = 0;
    for (auto& pixel1 : pixels) {
      const auto& neighbors47 = pixel1.neighbors_47;

      /*
      bool found = false;
      for (std::size_t i = 0; i < neighbors47.size(); ++i) {
        const std::size_t neighbor = neighbors47[i];
        if (neighbor == std::size_t(-1)) continue;
        const auto& pixel2 = pixels[neighbor];
        if (pixel2.binary == std::size_t(-1)) continue;
        if (pixel2.binary != pixel1.binary)
          found = true;
      }
      if (!found) continue; */

      for (std::size_t i = 0; i < neighbors47.size(); ++i) {
        const std::size_t neighbor = neighbors47[i];
        if (neighbor == std::size_t(-1)) continue;
        const auto& pixel2 = pixels[neighbor];
        if (pixel2.binary == std::size_t(-1)) continue;

        const auto& p = pixel1.point;
        const auto& q = pixel2.point;
        const auto  m = internal::middle_point_2(p, q);

        ++count;
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

    /* std::cout << "sizes: " << count << " " << dual.size() << std::endl; */
    /*
    save_dual_grid(
      dual, Color(0, 125, 0), lp_index); */
  }

  void apply_contouring(const std::size_t lp_index = 0) {

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
    std::vector<Segment_2> contour;
    contour.reserve(segments.size());
    for (const auto& segment : segments)
      contour.push_back(Segment_2(segment.source(), segment.target()));

    Saver saver;
    saver.save_polylines(
      contour, "/Users/monet/Documents/gf/lod/logs/buildings/tmp/duals/ms_contour-" +
      std::to_string(lp_index)); */
  }

  std::size_t get_cell_idx(const Indices& ns) {
    CGAL_assertion(ns.size() == 4);

    const std::size_t i0 = ns[0];
    const std::size_t i1 = ns[1];
    const std::size_t i2 = ns[2];
    const std::size_t i3 = ns[3];

    /* std::cout << i0 << " " << i1 << " " << i2 << " " << i3 << std::endl; */

    CGAL_assertion(i0 != std::size_t(-1));
    CGAL_assertion(i1 != std::size_t(-1));
    CGAL_assertion(i2 != std::size_t(-1));
    CGAL_assertion(i3 != std::size_t(-1));

    const std::size_t b0 = pixels[i0].binary;
    const std::size_t b1 = pixels[i1].binary;
    const std::size_t b2 = pixels[i2].binary;
    const std::size_t b3 = pixels[i3].binary;

    if (
      b0 == std::size_t(-1) ||
      b1 == std::size_t(-1) ||
      b2 == std::size_t(-1) ||
      b3 == std::size_t(-1) ) {

      return std::size_t(-1);
    }

    /* std::cout << b0 << " " << b1 << " " << b2 << " " << b3 << std::endl; */

    CGAL_assertion(b0 != std::size_t(-1));
    CGAL_assertion(b1 != std::size_t(-1));
    CGAL_assertion(b2 != std::size_t(-1));
    CGAL_assertion(b3 != std::size_t(-1));

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
    segments.push_back(segment);
  }

  void add_segment(
    const Pixel& px0, const Pixel& px1, const Pixel& px2, const Pixel& px3) {

    My_segment segment;
    segment.source_ = internal::middle_point_2(px0.point, px1.point);
    segment.target_ = internal::middle_point_2(px2.point, px3.point);
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

    Segment_wrappers wrappers;
    create_segment_wrappers(wrappers);

    using Region_growing = internal::Region_growing<
    Segment_wrappers, Oriented_neighbor_query, Oriented_image_region, Seed_map>;

    Indices seeds, idx_map;
    seeds.resize(wrappers.size());
    idx_map.resize(wrappers.size());

    for (std::size_t i = 0; i < wrappers.size(); ++i) {
      seeds[i] = wrappers[i].index;
      idx_map[wrappers[i].index] = i;
    }

    Seed_map seed_map(seeds);
    Oriented_neighbor_query onq(wrappers, idx_map);
    Oriented_image_region oir(wrappers, idx_map, 1);
    Region_growing region_growing(
      wrappers, onq, oir, seed_map);

    std::vector<Indices> regions;
    region_growing.detect(std::back_inserter(regions));

    /* std::cout << "num contours: " << regions.size() << std::endl; */

    Contour contour;
    contours.clear();

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
    Segment_wrappers& wrappers) {

    wrappers.clear();
    wrappers.resize(segments.size());

    Indices ns, nt;
    for (std::size_t i = 0; i < segments.size(); ++i) {
      wrappers[i].index = i;
      wrappers[i].neighbors.clear();

      const auto& source = segments[i].source();
      const auto& target = segments[i].target();
      find_neighbors(i, source, ns);
      find_neighbors(i, target, nt);

      add_neighbors(ns, wrappers[i].neighbors);
      add_neighbors(nt, wrappers[i].neighbors);
    }

    std::sort(wrappers.begin(), wrappers.end(),
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
    if (region.size() < 2) {

      const auto& curr = segments[region[0]];
      mp.point_ = curr.source();
      contour.points.push_back(mp);
      mp.point_ = curr.target();
      contour.points.push_back(mp);
      return;
    }

    const std::size_t rs = region.size() - 1;
    for (std::size_t i = 0; i < rs; ++i) {
      const std::size_t ip = i + 1;

      const auto& curr = segments[region[i]];
      const auto& next = segments[region[ip]];

      if (
        internal::are_equal_points_2(curr.source(), next.source()) ||
        internal::are_equal_points_2(curr.source(), next.target()) ) {

        mp.point_ = curr.target();
        contour.points.push_back(mp); continue;
      }

      if (
        internal::are_equal_points_2(curr.target(), next.source()) ||
        internal::are_equal_points_2(curr.target(), next.target()) ) {

        mp.point_ = curr.source();
        contour.points.push_back(mp); continue;
      }
    }

    const auto& curr = segments[region[rs - 1]];
    const auto& next = segments[region[rs]];

    if (internal::are_equal_points_2(curr.source(), next.source()) ) {

      mp.point_ = curr.source();
      contour.points.push_back(mp);
      mp.point_ = next.target();
      contour.points.push_back(mp);
      return;
    }

    if (internal::are_equal_points_2(curr.source(), next.target()) ) {

      mp.point_ = curr.source();
      contour.points.push_back(mp);
      mp.point_ = next.source();
      contour.points.push_back(mp);
      return;
    }

    if (internal::are_equal_points_2(curr.target(), next.source()) ) {

      mp.point_ = curr.target();
      contour.points.push_back(mp);
      mp.point_ = next.target();
      contour.points.push_back(mp);
      return;
    }

    if (internal::are_equal_points_2(curr.target(), next.target()) ) {

      mp.point_ = curr.target();
      contour.points.push_back(mp);
      mp.point_ = next.source();
      contour.points.push_back(mp);
      return;
    }
  }

  void type_contours() {

    for (auto& contour : contours) {
      if (line_found)
        set_linear_points(contour);
      type_contour(contour);
    }
  }

  void set_linear_points(Contour& contour) {

    auto& items = contour.points;
    for (auto& item : items) {

      const auto& p = item.point();
      const auto  q = line.projection(p);

      const FT distance = internal::distance(p, q);
      CGAL_assertion(noise_level_2 != FT(-1));
      if (distance < noise_level_2)
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

  void type_contour(Contour& contour) {

    Indices free, linear;
    auto& items = contour.points;
    for (std::size_t i = 0; i < items.size(); ++i) {

      // Other points.
      if (!items[i].belongs_to_line) {
        // Handle linear.
        if (linear.size() != 0)
          type_linear(linear, items);
        free.push_back(i);
      } else {
        // Handle free.
        if (free.size() != 0)
          type_free(free, items);
        linear.push_back(i);
      }
    }

    if (free.size() != 0)
      type_free(free, items);
    if (linear.size() != 0)
      type_linear(linear, items);

    // Fix unqiue linear points.
    for (std::size_t i = 0; i < items.size(); ++i) {
      if (items[i].int_type == Point_type::UNIQUE_LINEAR) {
        items[i].int_type = Point_type::LINEAR;
        continue;
      }
    }

    // Fix unique free points.
    for (std::size_t i = 0; i < items.size(); ++i) {
      if (items[i].int_type == Point_type::UNIQUE_FREE) {
        if (i > 1) {
          if (items[i - 1].int_type == Point_type::LINEAR) {
            items[i].int_type = Point_type::LINEAR;
            continue;
          }
        }
        if (i < items.size() - 1) {
          if (items[i + 1].int_type == Point_type::LINEAR) {
            items[i].int_type = Point_type::LINEAR;
            continue;
          }
        }
      }
    }
  }

  void type_linear(
    Indices& linear,
    std::vector<My_point>& items) {

    if (linear.size() >= 2) {
      for (const std::size_t idx : linear)
        items[idx].int_type = Point_type::LINEAR;
    } else {
      if (linear.size() == 1)
        items[linear[0]].int_type = Point_type::UNIQUE_LINEAR;
    }
    linear.clear();
  }

  void type_free(
    Indices& free,
    std::vector<My_point>& items) {

    if (free.size() >= 2) {
      for (const std::size_t idx : free)
        items[idx].int_type = Point_type::FREE;
    } else {
      if (free.size() == 1)
        items[free[0]].int_type = Point_type::UNIQUE_FREE;
    }
    free.clear();
  }

  void set_labels() {

    CGAL_assertion(label_pairs.size() == 1);
    const auto& label_pair = label_pairs[0];
    for (auto& contour : contours) {
      auto& items = contour.points;
      for (auto& item : items) {
        item.labels.clear();
        item.labels.insert(label_pair.first);
        item.labels.insert(label_pair.second);
      }
    }
  }

  void set_labels(
    const Size_pair& label_pair) {

    for (auto& contour : contours) {
      auto& items = contour.points;
      for (auto& item : items) {
        item.labels.clear();
        item.labels.insert(label_pair.first);
        item.labels.insert(label_pair.second);
      }
    }
  }

  void save_original_grid(
    const std::vector<Pixel>& pxs,
    const Color color0,
    const Color color1,
    const std::size_t lp_index = 0) {

    std::vector<Point_3> points0, points1;
    points0.reserve(pxs.size());
    points1.reserve(pxs.size());

    for (const auto& px : pxs) {
      const auto& point = px.point;

      if (px.binary == 0)
        points0.push_back(Point_3(point.x(), point.y(), FT(0)));
      else if (px.binary == 1)
        points1.push_back(Point_3(point.x(), point.y(), FT(0)));
    }
    Saver saver;
    saver.export_points(points0, color0,
    "/Users/monet/Documents/gf/lod/logs/buildings/tmp/duals/ms_original0-" +
    std::to_string(lp_index));
    saver.export_points(points1, color1,
    "/Users/monet/Documents/gf/lod/logs/buildings/tmp/duals/ms_original1-" +
    std::to_string(lp_index));
  }

  void save_dual_grid(
    const std::vector<Pixel>& pxs,
    const Color color,
    const std::size_t lp_index) {

    std::vector<Point_3> points;
    points.reserve(pxs.size());
    for (const auto& px : pxs) {
      const auto& point = px.point;
      points.push_back(Point_3(point.x(), point.y(), FT(0)));
    }
    Saver saver;
    saver.export_points(points, color,
    "/Users/monet/Documents/gf/lod/logs/buildings/tmp/duals/ms_dual-" +
    std::to_string(lp_index));
  }
};

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_INTERNAL_IMAGE_H
