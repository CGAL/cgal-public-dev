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

// Internal includes.
#include <CGAL/Levels_of_detail/internal/utils.h>
#include <CGAL/Levels_of_detail/internal/struct.h>

// Other includes.
#include <CGAL/Levels_of_detail/internal/Shape_detection/Region_growing.h>
#include <CGAL/Levels_of_detail/internal/Shape_detection/Linear_image_region.h>
#include <CGAL/Levels_of_detail/internal/Shape_detection/Planar_image_region.h>

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
  using Partition_edge_2 = internal::Partition_edge_2<Traits>;
  using Partition_face_2 = internal::Partition_face_2<Traits>;

  using Partition_2   = internal::Partition_2<Traits>;
  using LF_circulator = typename Triangulation::Delaunay::Line_face_circulator;
  using Size_pair     = std::pair<std::size_t, std::size_t>;
  using Idx_map       = std::map<Size_pair, std::size_t>;
  using Indices       = std::vector<std::size_t>;

  using Saver = Saver<Traits>;

  struct Pixel {
    Point_2 point;
    Indices neighbors_03;
    Indices neighbors_47;
    bool is_outer = false;
    std::size_t index = std::size_t(-1);
    std::size_t label = std::size_t(-1);
    std::size_t i = std::size_t(-1), j = std::size_t(-1);
  };

  struct Image {
    std::vector<Pixel> pixels;
    Indices seeds;
    std::vector<Size_pair> label_pairs;
    std::size_t num_labels;
    std::vector<Point_2> dual;
    bool is_ridge = false;
    
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
      neighbors = ns03;

      if (m_use_version_8) {
        const auto& ns47 = pixels[query_index].neighbors_47;
        for (const std::size_t neighbor : ns47)
          neighbors.push_back(neighbor);
      }
    }

    void clear() {
      pixels.clear();
      seeds.clear();
      label_pairs.clear();
      dual.clear();
      num_labels = 0;
    }

    void create_dual() {
      if (!is_ridge) return;
      std::vector<Segment_2> segments;
      create_contour(segments);
      make_dual(segments);
    }

  private:
    bool m_use_version_8 = false;

    void create_contour(
      std::vector<Segment_2>& segments) {
      segments.clear();

    }

    void make_dual(
      const std::vector<Segment_2>& segments) {
      dual.clear();

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
    const FT min_length_2,
    const FT angle_bound_2,
    const FT ordinate_bound_2) :
  m_boundary(boundary),
  m_lod0(lod0),
  m_image_ptr(image_ptr),
  m_partition_2(partition_2),
  m_min_length_2(min_length_2),
  m_angle_bound_2(angle_bound_2),
  m_ordinate_bound_2(ordinate_bound_2),
  m_pi(static_cast<FT>(CGAL_PI)) { 

    m_partition_2.clear();
    create_image();
  }

  void build() {

    std::size_t iter = 0;
    do {
      clean_image(); ++iter;
    } while (iter != 2);
    create_label_pairs();
    create_ridges();
    for (auto& ridge : m_ridges)
      ridge.create_dual();
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
  const FT m_min_length_2;
  const FT m_angle_bound_2;
  const FT m_ordinate_bound_2;
  const FT m_pi;

  Saver m_saver;
  Image m_image;
  std::vector<Image> m_ridges;

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
      "/Users/monet/Documents/lod/logs/buildings/tmp/image-clean.jpg", original);
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
    CGAL_assertion(best_idx != std::size_t(-1));
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
      "/Users/monet/Documents/lod/logs/buildings/tmp/ridges/ridge-" 
      + std::to_string(label_pair.first)  + "-" 
      + std::to_string(label_pair.second) + "-"
      + std::to_string(ridge_index) + ".jpg", image); 
  }
};

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_INTERNAL_PARTITION_BUILDER_FROM_IMAGE_2_H
