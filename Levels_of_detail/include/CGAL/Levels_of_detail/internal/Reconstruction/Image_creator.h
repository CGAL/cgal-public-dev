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

#ifndef CGAL_LEVELS_OF_DETAIL_INTERNAL_IMAGE_CREATOR_H
#define CGAL_LEVELS_OF_DETAIL_INTERNAL_IMAGE_CREATOR_H

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
#include <CGAL/point_generators_2.h>
#include <CGAL/Barycentric_coordinates_2/Segment_coordinates_2.h>

// Internal includes.
#include <CGAL/Levels_of_detail/internal/utils.h>
#include <CGAL/Levels_of_detail/internal/struct.h>

// Other includes.
#include <CGAL/Levels_of_detail/internal/Reconstruction/Image.h>
#include <CGAL/Levels_of_detail/internal/Shape_detection/Region_growing.h>
#include <CGAL/Levels_of_detail/internal/Shape_detection/Linear_image_region.h>
#include <CGAL/Levels_of_detail/internal/Shape_detection/Planar_image_region.h>

// Testing.
#include "../../../../../test/Levels_of_detail/include/Saver.h"

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

template<
typename GeomTraits,
typename ImagePointer>
class Image_creator {

public:
  using Traits = GeomTraits;
  using Image_ptr = ImagePointer;

  using FT = typename Traits::FT;
  using Point_2 = typename Traits::Point_2;
  using Segment_2 = typename Traits::Segment_2;
  using Line_2 = typename Traits::Line_2;
  using Intersect_2 = typename Traits::Intersect_2;

  using Image = internal::Image<Traits>;
  using Pixel = typename Image::Pixel;
  using Pixels = typename Image::Pixels;
  using My_point = typename Image::My_point;
  using Point_type = typename Image::Point_type;
  using Contour = typename Image::Contour;

  using Images = std::vector<Image>;

  using Indices = std::vector<std::size_t>;
  using Size_pair = std::pair<std::size_t, std::size_t>;
  using Idx_map = std::map<Size_pair, std::size_t>;

  using Seed_map = internal::Seed_property_map;
  using Linear_image_region = internal::Linear_image_region<Traits, Pixel>;
  using Planar_image_region = internal::Planar_image_region<Traits, Pixel>;

  using Saver = Saver<Traits>;

  using Point_pair = std::pair<Point_2, std::size_t>;
  using Point_pairs = std::vector<Point_pair>;
  using PS_generator = CGAL::Points_on_segment_2<Point_2>;
  using Point_pair_map = CGAL::First_of_pair_property_map<Point_pair>;
  using KNQ_pair =
    internal::K_neighbor_query<Traits, Point_pairs, Point_pair_map>;

  Image_creator(
    ImagePointer& image_ptr,
    const std::vector<Segment_2>& boundary,
    const FT noise_level_2) :
  m_image_ptr(image_ptr),
  m_boundary(boundary),
  m_noise_level_2(noise_level_2),
  m_pi(static_cast<FT>(CGAL_PI)) {

    create_boundary_knq();
  }

  void create_image() {
    build_image();
  }

  void clean_image() {
    std::size_t iter = 0;
    do { detect_and_clean(); ++iter;
    } while (iter != 2);
  }

  void create_label_pairs() {
    find_label_pairs();
  }

  void create_ridges_with_contours_v1() {
    create_ridges();
    create_contours();
  }

  void create_ridges_with_contours_v2() {

    m_ridges.clear();
    const auto& label_pairs = m_image.label_pairs;
    for (std::size_t i = 0; i < label_pairs.size(); ++i) {
      const auto& label_pair = label_pairs[i];
      create_ridges_from_label_pair_v1(i, label_pair);
    }
    clean_contours();
    save_ridge_polylines();
  }

  void create_ridges_with_contours_v3() {

    m_ridges.clear();
    const auto& label_pairs = m_image.label_pairs;
    for (std::size_t i = 0; i < label_pairs.size(); ++i) {
      const auto& label_pair = label_pairs[i];
      create_ridges_from_label_pair_v2(i, label_pair);
    }
    clean_contours();
    save_ridge_polylines();
  }

  void create_ridges() {

    m_ridges.clear();
    const auto& label_pairs = m_image.label_pairs;
    for (const auto& label_pair : label_pairs)
      extract_ridges(label_pair);

    /* std::cout << "num ridges: " << m_ridges.size() << std::endl; */
  }

  void create_contours() {

    for (auto& ridge : m_ridges) {
      ridge.create_line(
        m_image_ptr->get_plane_map());
      ridge.create_contours();
    }
    clean_contours();
    save_ridge_polylines();
  }

  const Images& get_ridges() const {
    return m_ridges;
  }

  const Image& get_image() const {
    return m_image;
  }

private:
  Image_ptr& m_image_ptr;
  const std::vector<Segment_2>& m_boundary;
  const FT m_noise_level_2;
  const FT m_pi;

  Image m_image;
  Images m_ridges;
  Point_pairs m_boundary_queries;
  std::shared_ptr<KNQ_pair> m_knq_ptr;

  void create_ridges_from_label_pair_v1(
    const std::size_t lp_index,
    const Size_pair& label_pair) {

    /* if (lp_index != 9) return; */
    Image ridge;
    ridge.clear();

    ridge.num_labels = m_image.num_labels;
    ridge.label_pairs.push_back(label_pair);
    ridge.is_ridge = true;
    ridge.noise_level_2 = m_noise_level_2;
    ridge.create_line(m_image_ptr->get_plane_map());

    m_image.noise_level_2 = m_noise_level_2;
    m_image.line = ridge.line;
    m_image.line_found = ridge.line_found;
    m_image.create_contours_from_label_pair_v1(
      lp_index, label_pair, ridge.contours);

    /* save_contours(ridge.contours, lp_index); */
    m_ridges.push_back(ridge);
  }

  void create_ridges_from_label_pair_v2(
    const std::size_t lp_index,
    const Size_pair& label_pair) {

    Image ridge;
    ridge.clear();

    ridge.num_labels = m_image.num_labels;
    ridge.label_pairs.push_back(label_pair);
    ridge.is_ridge = true;
    ridge.noise_level_2 = m_noise_level_2;
    ridge.create_line(m_image_ptr->get_plane_map());

    m_image.noise_level_2 = m_noise_level_2;
    m_image.line = ridge.line;
    m_image.line_found = ridge.line_found;
    m_image.create_contours_from_label_pair_v2(
      lp_index, label_pair, ridge.contours);

    /* save_contours(ridge.contours, lp_index); */
    m_ridges.push_back(ridge);
  }

  void save_contours(
    const std::vector<Contour>& contours,
    const std::size_t lp_index) {

    std::vector<Segment_2> segments;
    for (std::size_t i = 0; i < contours.size(); ++i) {
      const auto& contour = contours[i];
      if (contour.is_degenerated) continue;

      segments.clear();
      const auto& items = contour.points;
      for (std::size_t k = 0; k < items.size() - 1; ++k) {
        const std::size_t kp = k + 1;
        const auto& p = items[k].point();
        const auto& q = items[kp].point();
        segments.push_back(Segment_2(p, q));
      }

      Saver saver;
      saver.save_polylines(
        segments,
        "/Users/monet/Documents/gf/lod/logs/buildings/tmp/contours/ridge-contours-" +
        std::to_string(lp_index) + "-" + std::to_string(i));
    }
  }

  void create_boundary_knq() {

    create_boundary_queries();
    Point_pair_map pmap;
    m_knq_ptr =
      std::make_shared<KNQ_pair>(m_boundary_queries, FT(1), pmap);
  }

  void create_boundary_queries() {

    m_boundary_queries.clear();
    std::vector<Point_2> samples;
    const std::size_t num_samples_per_segment = 40;

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

  void build_image() {

    m_image.clear();
    const auto& original = m_image_ptr->get_image();
    const std::size_t num_labels = m_image_ptr->get_num_labels();
    m_image.num_labels = num_labels;
    auto& pixels = m_image.pixels;
    auto& seeds = m_image.seeds;

    Pixel ipixel;
    Idx_map idx_map;

    std::size_t index = 0;
    for (std::size_t i = 0; i < original.rows; ++i) {
      for (std::size_t j = 0; j < original.cols; ++j) {

        const auto& cell = original.grid[i][j];
        const std::size_t label = m_image_ptr->get_label(
          cell.zr, cell.zg, cell.zb);

        if (label == num_labels) {
          ipixel.label = std::size_t(-1);
          ipixel.is_outer = true;
        } else {
          ipixel.label = label;
          ipixel.is_outer = false;
        }

        ipixel.i = i; ipixel.j = j;
        ipixel.index = index;
        ipixel.point = m_image_ptr->get_point(i, j);
        pixels.push_back(ipixel);

        if (ipixel.label == std::size_t(-1))
          seeds.push_back(ipixel.label);
        else
          seeds.push_back(ipixel.index);

        idx_map[std::make_pair(ipixel.i, ipixel.j)] = ipixel.index;
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

    /*
    const auto& pixel1 = pixel;
    for (const std::size_t idx : pixel.neighbors_03) {
      const auto& pixel2 = m_image.pixels[idx];

      if (pixel1.point == pixel2.point) {
        std::cout.precision(30);
        std::cout << "Image error: " << pixel1.neighbors_03.size() << std::endl;
        std::cout << pixel1.i << " " << pixel1.j << std::endl;
        std::cout << m_image_ptr->get_point(pixel1.i, pixel1.j) << std::endl;
        std::cout << pixel2.i << " " << pixel2.j << std::endl;
        std::cout << m_image_ptr->get_point(pixel2.i, pixel2.j) << std::endl;
        exit(EXIT_FAILURE);
      }
    } */
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

  void detect_and_clean() {

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

    /* std::cout << "num labeled regions: " << regions.size() << std::endl; */

    // find best label
    auto& original = m_image_ptr->get_image();

    m_image_ptr->save_image(
      "/Users/monet/Documents/gf/lod/logs/buildings/tmp/image-clean-0.jpg", original);

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
      "/Users/monet/Documents/gf/lod/logs/buildings/tmp/image-clean-1.jpg", original);

    // remove degenerated corner pixels
    for (auto& pixel : pixels) {
      const auto& neighbors = pixel.neighbors_03;

      const std::size_t result = pixels[neighbors[0]].label;
      if (result == std::size_t(-1))
        continue;

      std::size_t count = 1;
      for (std::size_t k = 1; k < neighbors.size(); ++k)
        if (pixels[neighbors[k]].label == result) ++count;

      if (count == neighbors.size()) {
        pixel.label = result;

        const std::size_t i = pixel.i;
        const std::size_t j = pixel.j;

        auto& cell = original.grid[i][j];
        const auto& p = m_image_ptr->get_label_map().at(pixel.label);
        cell.zr = p.x(); cell.zg = p.y(); cell.zb = p.z();
      }
    }

    for (auto& pixel : pixels) {
      if (pixel.label != std::size_t(-1)) continue;
      const auto& neighbors = pixel.neighbors_03;

      std::size_t count = 0;
      for (std::size_t k = 0; k < neighbors.size(); ++k)
        if (pixels[neighbors[k]].label != std::size_t(-1)) ++count;

      if (count == neighbors.size()) {
        pixel.label = pixels[neighbors[0]].label;

        const std::size_t i = pixel.i;
        const std::size_t j = pixel.j;

        auto& cell = original.grid[i][j];
        const auto& p = m_image_ptr->get_label_map().at(pixel.label);
        cell.zr = p.x(); cell.zg = p.y(); cell.zb = p.z();
      }
    }

    for (auto& pixel : pixels)
      pixel.used = false;

    m_image_ptr->save_image(
      "/Users/monet/Documents/gf/lod/logs/buildings/tmp/image-clean-2.jpg", original);

    // remove degenerated internal pixels
    for (auto& pixel : pixels) {
      if (pixel.used) continue;

      pixel.used = true;
      const auto& neighbors03 = pixel.neighbors_03;
      const auto& neighbors47 = pixel.neighbors_47;

      std::size_t result = std::size_t(-1);
      const bool is_degenerated =
        is_degenerated_corner(pixel, neighbors03, neighbors47, result);

      if (is_degenerated) {
        pixel.label = result;

        const std::size_t i = pixel.i;
        const std::size_t j = pixel.j;

        auto& cell = original.grid[i][j];
        const auto& p = m_image_ptr->get_label_map().at(pixel.label);
        cell.zr = p.x(); cell.zg = p.y(); cell.zb = p.z();
      }
    }

    for (auto& pixel : pixels)
      pixel.used = false;

    m_image_ptr->save_image(
      "/Users/monet/Documents/gf/lod/logs/buildings/tmp/image-clean.jpg", original);
  }

  bool is_degenerated_corner(
    const Pixel& query,
    const Indices& neighbors03,
    const Indices& neighbors47,
    std::size_t& result) {

    result = std::size_t(-1);
    if (query.label == std::size_t(-1))
      return false;

    const auto& pixels = m_image.pixels;
    if (is_case_0(query, neighbors03, neighbors47)) {
      const std::size_t l1 = pixels[neighbors03[3]].label;
      const std::size_t l2 = pixels[neighbors03[0]].label;
      if (l1 != std::size_t(-1)) result = l1;
      else if (l2 != std::size_t(-1)) result = l2;
      else return false;
      return true;
    }
    if (is_case_1(query, neighbors03, neighbors47)) {
      const std::size_t l1 = pixels[neighbors03[0]].label;
      const std::size_t l2 = pixels[neighbors03[1]].label;
      if (l1 != std::size_t(-1)) result = l1;
      else if (l2 != std::size_t(-1)) result = l2;
      else return false;
      return true;
    }
    if (is_case_2(query, neighbors03, neighbors47)) {
      const std::size_t l1 = pixels[neighbors03[1]].label;
      const std::size_t l2 = pixels[neighbors03[2]].label;
      if (l1 != std::size_t(-1)) result = l1;
      else if (l2 != std::size_t(-1)) result = l2;
      else return false;
      return true;
    }
    if (is_case_3(query, neighbors03, neighbors47)) {
      const std::size_t l1 = pixels[neighbors03[2]].label;
      const std::size_t l2 = pixels[neighbors03[3]].label;
      if (l1 != std::size_t(-1)) result = l1;
      else if (l2 != std::size_t(-1)) result = l2;
      else return false;
      return true;
    }
    return false;
  }

  bool is_case_0(
    const Pixel& query,
    const Indices& neighbors03,
    const Indices& neighbors47) {

    auto& pixels = m_image.pixels;

    const std::size_t ref = query.label;
    const std::size_t other = pixels[neighbors47[0]].label;
    if (ref != other) return false;

    const std::size_t l1 = pixels[neighbors03[3]].label;
    const std::size_t l2 = pixels[neighbors03[0]].label;
    if (ref == l1 || ref == l2) return false;

    pixels[neighbors03[3]].used = true;
    pixels[neighbors03[0]].used = true;
    pixels[neighbors47[0]].used = true;
    return true;
  }

  bool is_case_1(
    const Pixel& query,
    const Indices& neighbors03,
    const Indices& neighbors47) {

    auto& pixels = m_image.pixels;

    const std::size_t ref = query.label;
    const std::size_t other = pixels[neighbors47[1]].label;
    if (ref != other) return false;

    const std::size_t l1 = pixels[neighbors03[0]].label;
    const std::size_t l2 = pixels[neighbors03[1]].label;
    if (ref == l1 || ref == l2) return false;

    pixels[neighbors03[0]].used = true;
    pixels[neighbors03[1]].used = true;
    pixels[neighbors47[1]].used = true;
    return true;
  }

  bool is_case_2(
    const Pixel& query,
    const Indices& neighbors03,
    const Indices& neighbors47) {

    auto& pixels = m_image.pixels;

    const std::size_t ref = query.label;
    const std::size_t other = pixels[neighbors47[2]].label;
    if (ref != other) return false;

    const std::size_t l1 = pixels[neighbors03[1]].label;
    const std::size_t l2 = pixels[neighbors03[2]].label;
    if (ref == l1 || ref == l2) return false;

    pixels[neighbors03[1]].used = true;
    pixels[neighbors03[2]].used = true;
    pixels[neighbors47[3]].used = true;
    return true;
  }

  bool is_case_3(
    const Pixel& query,
    const Indices& neighbors03,
    const Indices& neighbors47) {

    auto& pixels = m_image.pixels;

    const std::size_t ref = query.label;
    const std::size_t other = pixels[neighbors47[3]].label;
    if (ref != other) return false;

    const std::size_t l1 = pixels[neighbors03[2]].label;
    const std::size_t l2 = pixels[neighbors03[3]].label;
    if (ref == l1 || ref == l2) return false;

    pixels[neighbors03[2]].used = true;
    pixels[neighbors03[3]].used = true;
    pixels[neighbors47[3]].used = true;
    return true;
  }

  std::size_t get_best_label(
    const Indices& region) {

    const auto& pixels = m_image.pixels;
    const std::size_t num_labels = m_image.num_labels;

    Indices neighbors;
    Indices nums(num_labels, 0);
    for (const std::size_t idx : region) {

      const std::size_t ref_label = pixels[idx].label;
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
    return best_idx;
  }

  void find_label_pairs() {

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
        if (pixels[neighbor].label != pixel.label) {

          const std::size_t l1 = pixel.label;
          const std::size_t l2 = pixels[neighbor].label;

          if (l2 == std::size_t(-1)) {
            unique.insert(std::make_pair(l1, l2));
            continue;
          }

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

    /* std::cout << "num label pairs: " << label_pairs.size() << std::endl; */
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
      Pixels, Image, Linear_image_region>;

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
    ridge.clear();

    ridge.num_labels = num_labels;
    ridge.label_pairs = label_pairs;
    ridge.is_ridge = true;
    ridge.noise_level_2 = m_noise_level_2;

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
    m_ridges.push_back(ridge);

    /* save_ridge_image(ridge.label_pairs[0], rindex, ridge.pixels); */
  }

  void clean_contours() {

    set_types();

    /* // removed in version 2
    mark_internal_degenerated_contours_before(); */

    mark_end_points();

    /* // removed in version 1
    for (auto& ridge : m_ridges) {
      for (auto& contour : ridge.contours) {
        if (contour.skip()) continue;
        contour.truncate();
      }
    } */

    intersect_contours();
    update_internal_points();

    /* // removed in version 2
    mark_internal_degenerated_contours_after(); */

    /* // removed in version 3
    type_end_points(); */

    /* // removed in version 3
    update_boundary_points(); */

    /* // removed in version 2
    mark_boundary_degenerated_contours(); */

    /* // removed in version 3
    snap_contours(); */

    /* // removed in version 2
    remove_duplicated_edges(); */
  }

  void update_internal_points() {
    for (const auto& ridge : m_ridges) {
      for (const auto& contour : ridge.contours) {
        for (const auto& item : contour.points) {
          const auto& p = item.point_;
          update_closest_point(p);
        }
      }
    }
  }

  void update_closest_point(
    const Point_2& query) {

    for (auto& ridge : m_ridges) {
      for (auto& contour : ridge.contours) {
        for (auto& item : contour.points) {
          if (
            internal::are_equal_points_2(query, item.point_) &&
            query != item.point_) {
            item.point_ = query;
          }
        }
      }
    }
  }

  void set_types() {

    for (auto& ridge : m_ridges) {
      for (auto& contour : ridge.contours) {

        auto& items = contour.points;
        for (auto& item : items)
          item.end_type = item.int_type;
      }
    }
  }

  void mark_internal_degenerated_contours_before() {

    for (std::size_t i = 0; i < m_ridges.size(); ++i) {
      for (auto& contour : m_ridges[i].contours) {
        if (contour.is_closed) continue;
        if (contour.points.size() >= 10) continue;

        if (belongs_to_other_contours(i, contour.points))
          contour.is_degenerated = true;
      }
    }
  }

  bool belongs_to_other_contours(
    const std::size_t skip, const std::vector<My_point>& queries) {

    std::size_t count = 0;
    for (const auto& q : queries) {
      for (std::size_t i = 0; i < m_ridges.size(); ++i) {
        if (i == skip) continue;

        bool found_out = false;
        for (const auto& contour : m_ridges[i].contours) {
          if (contour.skip()) continue;

          bool found_in = false;
          for (const auto& p : contour.points) {
            if (internal::are_equal_points_2(p.point(), q.point())) {
              count++; found_in = true; break;
            }
          }
          if (found_in) {
            found_out = true; break;
          }
        }
        if (found_out) break;
      }
    }
    return ( count == queries.size() );
  }

  void mark_end_points() {

    /*
    FT max_dist = FT(0);
    for (const auto& ridge : m_ridges) {

      bool found = false;
      for (const auto& contour : ridge.contours) {
        if (contour.points.size() < 2) continue;

        const auto& p = contour.points[0];
        const auto& q = contour.points[1];
        max_dist = internal::distance(p.point(), q.point());
        found = true; break;
      }
      if (found) break;
    }
    max_dist *= FT(2); */
    const FT max_dist = internal::tolerance<FT>() * FT(2);

    for (std::size_t i = 0; i < m_ridges.size(); ++i) {
      for (auto& contour : m_ridges[i].contours) {
        if (contour.skip()) continue;
        for (auto& p : contour.points)
          add_contour_neighbors(max_dist, i, p);
      }
    }
  }

  void add_contour_neighbors(
    const FT max_dist,
    const std::size_t skip,
    My_point& p) {

    for (std::size_t i = 0; i < m_ridges.size(); ++i) {
      if (i == skip) continue;
      for (std::size_t j = 0; j < m_ridges[i].contours.size(); ++j) {
        const auto& contour = m_ridges[i].contours[j];
        if (contour.skip()) continue;

        const auto& q1 = contour.points[0];
        const auto& q2 = contour.points[contour.points.size() - 1];
        const FT dist1 = internal::distance(p.point(), q1.point());
        const FT dist2 = internal::distance(p.point(), q2.point());

        if (dist1 < dist2) {
          if (dist1 <= max_dist) {
            p.neighbors.insert(std::make_pair(i, j));
          }
        } else {
          if (dist2 <= max_dist) {
            p.neighbors.insert(std::make_pair(i, j));
          }
        }

        /*
        for (const auto& q : contour.points) {
          if (internal::are_equal_points_2(p.point(), q.point()))
            p.neighbors.insert(std::make_pair(i, j));
        } */
      }
    }
  }

  void intersect_contours() {

    const FT max_dist = internal::tolerance<FT>() * FT(2);
    for (auto& ridge : m_ridges) {
      for (auto& contour : ridge.contours) {
        if (contour.skip()) continue;

        auto& items = contour.points;
        const std::size_t nump = items.size();

        auto& p = items[0];
        auto& q = items[nump - 1];

        update_end_point(max_dist, p);
        update_end_point(max_dist, q);
      }
    }

    for (std::size_t i = 0; i < m_ridges.size(); ++i) {
      for (auto& contour : m_ridges[i].contours) {
        if (contour.skip()) continue;

        auto& items = contour.points;
        const std::size_t nump = items.size();

        auto& p = items[0];
        auto& q = items[nump - 1];

        average_point(max_dist, i, p);
        average_point(max_dist, i, q);
      }
    }
  }

  void update_end_point(
    const FT max_dist,
    My_point& query) {

    FT x = query.point().x();
    FT y = query.point().y();

    const auto& neighbors = query.neighbors;
    if (neighbors.size() == 0) return;

    for (const auto& pair : neighbors) {
      const auto& contour = m_ridges[pair.first].contours[pair.second];

      const auto& items = contour.points;
      const std::size_t nump = items.size();

      const auto& p = contour.points[0];
      const auto& q = contour.points[nump - 1];

      const FT dist1 = internal::distance(query.point(), p.point());
      const FT dist2 = internal::distance(query.point(), q.point());

      if (dist1 < dist2) {
        if (dist1 <= max_dist) {
          x += p.point().x(); y += p.point().y();
        }
      } else {
        if (dist2 <= max_dist) {
          x += q.point().x(); y += q.point().y();
        }
      }
    }
    x /= static_cast<FT>(neighbors.size() + 1);
    y /= static_cast<FT>(neighbors.size() + 1);

    const Point_2 center = Point_2(x, y);
    query.point_ = center;

    for (const auto& pair : neighbors) {
      auto& contour = m_ridges[pair.first].contours[pair.second];

      auto& items = contour.points;
      const std::size_t nump = items.size();

      auto& p = contour.points[0];
      auto& q = contour.points[nump - 1];

      const FT dist1 = internal::distance(query.point(), p.point());
      const FT dist2 = internal::distance(query.point(), q.point());

      if (dist1 < dist2) {
        if (dist1 <= max_dist) {
          p.point_ = center;
        }
      } else {
        if (dist2 <= max_dist) {
          q.point_ = center;
        }
      }
    }
  }

  void average_point(
    const FT max_dist,
    const std::size_t skip,
    My_point& query) {

    const auto& neighbors = query.neighbors;
    if (neighbors.size() == 0) return;

    for (std::size_t i = 0; i < m_ridges.size(); ++i) {
      if (i == skip) continue;

      for (auto& contour : m_ridges[i].contours) {
        if (contour.skip()) continue;

        auto& items = contour.points;
        const std::size_t nump = items.size();

        auto& p = items[0];
        auto& q = items[nump - 1];

        const FT dist1 = internal::distance(query.point(), p.point());
        const FT dist2 = internal::distance(query.point(), q.point());

        if (dist1 < dist2) {
          if (dist1 <= max_dist) {
            p.point_ = query.point_;
          }
        } else {
          if (dist2 <= max_dist) {
            q.point_ = query.point_;
          }
        }
      }
    }
  }

  void mark_internal_degenerated_contours_after() {
    for (auto& ridge : m_ridges) {
      for (auto& contour : ridge.contours) {
        if (contour.points.size() < 2) {
          contour.is_degenerated = true; continue;
        }

        FT dist = FT(0);
        for (std::size_t i = 0; i < contour.points.size() - 1; ++i) {
          const auto& p = contour.points[i];
          const auto& q = contour.points[i + 1];
          dist += internal::distance(p.point(), q.point());
        }

        if (dist <= internal::tolerance<FT>()) {
          std::cout << "degenerated contour removed" << std::endl;
          contour.is_degenerated = true;
        }
      }
    }
  }

  void type_end_points() {

    for (auto& ridge : m_ridges) {
      for (auto& contour : ridge.contours) {
        if (contour.skip()) continue;

        auto& items = contour.points;
        const std::size_t nump = items.size();

        auto& p = items[0];
        auto& q = items[nump - 1];

        if (p.neighbors.size() < 2)
          p.end_type = Point_type::BOUNDARY;
        else
          p.end_type = Point_type::CORNER;

        if (q.neighbors.size() < 2)
          q.end_type = Point_type::BOUNDARY;
        else
          q.end_type = Point_type::CORNER;
      }
    }
  }

  void update_boundary_points() {

    for (auto& ridge : m_ridges) {
      for (auto& contour : ridge.contours) {
        if (contour.skip()) continue;

        auto& items = contour.points;
        if (items[0].end_type == Point_type::BOUNDARY) {
          const bool success = intersect_with_boundary(items[0]);
          /* if (success) clean_intersection_begin(items); */
        }
        if (items[items.size() - 1].end_type == Point_type::BOUNDARY) {
          const bool success = intersect_with_boundary(items[items.size() - 1]);
          /* if (success) clean_intersection_end(items); */
        }
      }
    }
  }

  bool intersect_with_boundary(My_point& query) {

    Indices closest;
    (*m_knq_ptr)(query.point(), closest);
    const std::size_t bd_idx = m_boundary_queries[closest[0]].second;

    const auto& segment = m_boundary[bd_idx];
    const auto& s = segment.source();
    const auto& t = segment.target();
    const Line_2 line2 = Line_2(s, t);

    const auto proj = line2.projection(query.point());
    if (internal::distance(query.point(), proj) >= m_noise_level_2) {
      query.end_type = Point_type::DISTANT_BOUNDARY; return false;
    }

    query.bd_idx = bd_idx;
    query.point_ = proj;

    return true;
  }

  void clean_intersection_begin(
    std::vector<My_point>& items) {

    const std::size_t nump = items.size();
    const std::size_t nums = ( nump >= 6 ? 5 : nump - 1 );
    const auto& query = items[0];

    const std::size_t bd_idx = query.bd_idx;
    const auto& segment = m_boundary[bd_idx];
    const auto& s = segment.source();
    const auto& t = segment.target();
    const Line_2 ref_line = Line_2(s, t);

    Point_2 start_point;
    std::size_t start = std::size_t(-1);

    for (std::size_t i = 1; i < nums - 1; ++i) {
      const std::size_t ip = i + 1;

      const auto& p = items[i];
      const auto& q = items[ip];

      const Line_2 line = Line_2(p.point(), q.point());
      Point_2 r;
      const bool success = intersect_2(line, ref_line, r);
      if (success) {
        if (belongs_to_segment(s, r, t)) {
          start = ip; start_point = ref_line.projection(p.point()); continue;
        }
      }
    }

    if (start != std::size_t(-1)) {
      std::vector<My_point> clean;

      for (std::size_t k = start; k < items.size(); ++k)
        clean.push_back(items[k]);
      clean[0].end_type = Point_type::BOUNDARY;
      clean[0].point_ = start_point;
      clean[0].bd_idx = bd_idx;
      items.clear(); items = clean;
    }
  }

  void clean_intersection_end(
    std::vector<My_point>& items) {

    const std::size_t nump = items.size();
    const std::size_t nums = ( nump >= 6 ? 5 : nump - 1 );
    const auto& query = items[nump - 1];

    const std::size_t bd_idx = query.bd_idx;
    const auto& segment = m_boundary[bd_idx];
    const auto& s = segment.source();
    const auto& t = segment.target();
    const Line_2 ref_line = Line_2(s, t);

    Point_2 end_point;
    std::size_t end = std::size_t(-1);

    for (std::size_t i = nump - nums + 1; i < nump - 1; ++i) {
      const std::size_t ip = i + 1;

      const auto& p = items[i];
      const auto& q = items[ip];

      const Line_2 line = Line_2(p.point(), q.point());
      Point_2 r;
      const bool success = intersect_2(line, ref_line, r);
      if (success) {
        if (belongs_to_segment(s, r, t)) {
          end = ip; end_point = ref_line.projection(p.point()); break;
        }
      }
    }

    if (end != std::size_t(-1)) {
      std::vector<My_point> clean;

      for (std::size_t k = 0; k < end; ++k)
        clean.push_back(items[k]);
      const std::size_t numc = clean.size();
      clean[numc - 1].end_type = Point_type::BOUNDARY;
      clean[numc - 1].point_ = end_point;
      clean[numc - 1].bd_idx = bd_idx;
      items.clear(); items = clean;
    }
  }

  bool belongs_to_segment(
    const Point_2& s, const Point_2& r, const Point_2& t) {

    const auto res =
      CGAL::Barycentric_coordinates::compute_segment_coordinates_2(
        s, t, r, Traits());
    return ( res[0] >= FT(0) && res[1] >= FT(0) );
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

  void mark_boundary_degenerated_contours() {

    for (std::size_t i = 0; i < m_ridges.size(); ++i) {
      for (auto& contour : m_ridges[i].contours) {
        if (contour.skip()) continue;

        if (belongs_to_boundary(contour.points))
          contour.is_degenerated = true;
      }
    }
  }

  bool belongs_to_boundary(
    const std::vector<My_point>& items) {

    const std::size_t nump = items.size();
    const auto& p = items[0];
    const auto& q = items[nump - 1];

    if (p.end_type != Point_type::BOUNDARY)
      return false;
    if (q.end_type != Point_type::BOUNDARY)
      return false;

    if (p.bd_idx == std::size_t(-1))
      return false;
    if (q.bd_idx == std::size_t(-1))
      return false;

    if (p.bd_idx != q.bd_idx)
      return false;

    const auto& segment = m_boundary[p.bd_idx];
    const auto& s = segment.source();
    const auto& t = segment.target();

    for (std::size_t i = 1; i < nump - 1; ++i) {
      const auto& r = items[i].point();
      if (!belongs_to_segment(s, r, t)) return false;
    }
    return true;
  }

  void snap_contours() {

    for (std::size_t i = 0; i < m_ridges.size(); ++i) {
      for (auto& contour : m_ridges[i].contours) {
        if (contour.skip()) continue;

        auto& items = contour.points;
        const std::size_t nump = items.size();
        auto& p = items[0];
        auto& q = items[nump - 1];

        if (p.end_type == Point_type::BOUNDARY)
          snap_boundary_point(p);
        if (q.end_type == Point_type::BOUNDARY)
          snap_boundary_point(q);
      }
    }
  }

  void snap_boundary_point(My_point& query) {

    const std::size_t bd_idx = query.bd_idx;
    const auto& segment = m_boundary[bd_idx];
    const auto& s = segment.source();
    const auto& t = segment.target();
    const auto  m = internal::middle_point_2(s, t);

    const FT dist1 = internal::distance(query.point(), s);
    const FT dist2 = internal::distance(query.point(), t);
    const FT dist3 = internal::distance(query.point(), m);

    const FT eps = m_noise_level_2 / FT(4);
    if (dist1 <= eps) {
      query.point_ = s; return;
    }
    if (dist2 <= eps) {
      query.point_ = t; return;
    }
    if (dist3 <= eps) {
      query.point_ = m; return;
    }
  }

  void remove_duplicated_edges() {

    bool found = false;
    for (std::size_t i = 0; i < m_ridges.size(); ++i) {
      for (auto& contour : m_ridges[i].contours) {
        if (contour.skip()) continue;

        auto& items = contour.points;
        for (std::size_t j = 0; j < items.size() - 1; ++j) {
          const std::size_t jp = j + 1;
          if (is_edge_on_another_contour(i, i, j, items[j], items[jp])) {
            found = true; std::cout << "duplicate detected" << std::endl;
          }
        }

        std::vector<My_point> clean;
        for (std::size_t j = 0; j < items.size(); ++j) {
          if (!items[j].is_removed) clean.push_back(items[j]);
        }
        items = clean;
      }
    }

    if (found) {
      std::vector<My_point> clean;
      for (std::size_t i = 0; i < m_ridges.size(); ++i) {
        for (auto& contour : m_ridges[i].contours) {
          if (contour.skip()) continue;

          clean.clear();
          auto& items = contour.points;
          for (std::size_t j = 0; j < items.size(); ++j) {
            if (!items[j].is_removed) clean.push_back(items[j]);
          }
          items = clean;
        }
      }
    }
  }

  bool is_edge_on_another_contour(
    const std::size_t skip,
    const std::size_t ii, const std::size_t jj,
    My_point& a, My_point& b) {

    for (std::size_t i = 0; i < m_ridges.size(); ++i) {
      if (i == skip) continue;

      for (auto& contour : m_ridges[i].contours) {
        if (contour.skip()) continue;

        auto& items = contour.points;
        for (std::size_t j = 0; j < items.size() - 1; ++j) {
          const std::size_t jp = j + 1;

          auto& c = items[j];
          auto& d = items[jp];

          if (
            internal::are_equal_points_2(a.point(), c.point()) &&
            internal::are_equal_points_2(b.point(), d.point()) ) {

            a.merge(c); b.merge(d);
            a.merge(b); c.merge(d);
            b.is_removed = true;
            d.is_removed = true;

            a.add_neighbor(i, j);
            c.add_neighbor(ii, jj);

            return true;
          }

          if (
            internal::are_equal_points_2(a.point(), d.point()) &&
            internal::are_equal_points_2(b.point(), c.point()) ) {

            a.merge(d); b.merge(c);
            a.merge(b); d.merge(c);
            b.is_removed = true;
            c.is_removed = true;

            a.add_neighbor(i, j);
            d.add_neighbor(ii, jj);

            return true;
          }
        }
      }
    }
    return false;
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

  void save_ridge_polylines(
    const bool save_all = true) {

    if (save_all)
      save_ridge_polylines_all_at_once();
    else
      save_ridge_polylines_one_by_one();
  }

  void save_ridge_polylines_one_by_one() {

    std::vector<Segment_2> segments;
    for (std::size_t i = 0; i < m_ridges.size(); ++i) {
      for (std::size_t j = 0; j < m_ridges[i].contours.size(); ++j) {
        const auto& contour = m_ridges[i].contours[j];
        if (contour.is_degenerated) continue;

        segments.clear();
        const auto& items = contour.points;
        for (std::size_t k = 0; k < items.size() - 1; ++k) {
          const std::size_t kp = k + 1;
          const auto& p = items[k].point();
          const auto& q = items[kp].point();
          segments.push_back(Segment_2(p, q));
        }

        Saver saver;
        saver.save_polylines(
          segments,
          "/Users/monet/Documents/gf/lod/logs/buildings/tmp/ridges/ridge-contours-" +
          std::to_string(i) + "-" + std::to_string(j));
      }
    }
  }

  void save_ridge_polylines_all_at_once() {

    std::vector<Segment_2> segments;
    /* const auto& ridge = m_ridges[0]; */

    for (const auto& ridge : m_ridges) {
      for (const auto& contour : ridge.contours) {
        if (contour.is_degenerated) continue;

        const auto& items = contour.points;
        for (std::size_t k = 0; k < items.size() - 1; ++k) {
          const std::size_t kp = k + 1;
          const auto& p = items[k].point();
          const auto& q = items[kp].point();
          segments.push_back(Segment_2(p, q));
        }
      }
    }

    Saver saver;
    saver.save_polylines(
      segments, "/Users/monet/Documents/gf/lod/logs/buildings/tmp/ridge-contours");
  }
};

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_INTERNAL_IMAGE_CREATOR_H
