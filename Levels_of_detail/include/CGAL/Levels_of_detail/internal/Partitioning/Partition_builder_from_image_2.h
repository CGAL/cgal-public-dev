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
#include <vector>

// CGAL includes.
#include <CGAL/property_map.h>

// Internal includes.
#include <CGAL/Levels_of_detail/internal/utils.h>
#include <CGAL/Levels_of_detail/internal/struct.h>

#include <CGAL/Levels_of_detail/internal/Shape_detection/Region_growing.h>
#include <CGAL/Levels_of_detail/internal/Shape_detection/Planar_image_region.h>
#include <CGAL/Levels_of_detail/internal/Spatial_search/Image_neighbor_query.h>

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

template<
typename GeomTraits,
typename ImagePointer>
class Partition_builder_from_image_2 {

public:
  using Traits = GeomTraits;

  using FT = typename Traits::FT;
  using Point_2 = typename Traits::Point_2;
  using Segment_2 = typename Traits::Segment_2;

  using Triangulation = internal::Triangulation<Traits>;
  using Face_handle = typename Triangulation::Delaunay::Face_handle;
  using Vertex_handle = typename Triangulation::Delaunay::Vertex_handle;
  using Partition_edge_2 = internal::Partition_edge_2<Traits>;
  using Partition_face_2 = internal::Partition_face_2<Traits>;

  using Partition_2 = internal::Partition_2<Traits>;
  using LF_circulator = typename Triangulation::Delaunay::Line_face_circulator;
  using Size_pair = std::pair<std::size_t, std::size_t>;

  struct Pixel {
    Pixel(
      const std::size_t index_, 
      const std::size_t i_,
      const std::size_t j_,
      const bool is_boundary_,
      const std::size_t label_) :
    index(index_),
    i(i_), j(j_),
    is_boundary(is_boundary_),
    label(label_)
    { }

    std::size_t index = std::size_t(-1);
    std::size_t i = std::size_t(-1);
    std::size_t j = std::size_t(-1);
    bool is_boundary = false;
    std::size_t label = std::size_t(-1);
  };

  using Image = std::vector<Pixel>;
  using Idx_map = std::map<Size_pair, std::size_t>;
  using Seeds = std::vector<std::size_t>;
  using Seed_map = internal::Seed_property_map;

  using Image_neighbor_query = internal::Image_neighbor_query<Traits, Pixel>;
  using Planar_image_region = internal::Planar_image_region<Traits, Pixel>;
  using Region_growing = internal::Region_growing<
    Seeds, Image_neighbor_query, Planar_image_region, Seed_map>;

  Partition_builder_from_image_2(
    const std::vector<Segment_2>& boundary,
    ImagePointer& image_ptr,
    Partition_2& partition_2) :
  m_boundary(boundary),
  m_image_ptr(image_ptr),
  m_partition_2(partition_2)
  { }

  void build() {
    m_partition_2.clear();

    m_constraints.clear();
    add_boundary_constraints();
    clean_image();
    add_internal_constraints();
    create_triangulation();
  }

private:
  const std::vector<Segment_2>& m_boundary;

  ImagePointer& m_image_ptr;
  Partition_2& m_partition_2;

  std::vector<Segment_2> m_constraints;
  Triangulation m_base;

  void add_boundary_constraints() {
    for (const auto& segment : m_boundary)
      m_constraints.push_back(segment);
  }

  void clean_image() {

    auto& original = m_image_ptr->get_image();
    
    Image image;
    Idx_map idx_map;
    Seeds seeds;

    std::size_t count = 0;
    for (std::size_t i = 0; i < original.rows; ++i) {
      for (std::size_t j = 0; j < original.cols; ++j) {

        if (
          i == 0 || j == 0 || 
          i == original.rows - 1 || j == original.cols - 1) {
          
          image.push_back(Pixel(count, i, j, true, std::size_t(-1)));
          seeds.push_back(std::size_t(-1));

        } else {

          const auto& cell = original.grid[i][j];
          const std::size_t label = m_image_ptr->get_label(
            cell.zr, cell.zg, cell.zb);
          
          if (label != m_image_ptr->get_num_labels()) {
            
            image.push_back(Pixel(count, i, j, false, label));
            seeds.push_back(count);

          } else {

            image.push_back(Pixel(count, i, j, true, std::size_t(-1)));
            seeds.push_back(std::size_t(-1));

          }
        }
        idx_map[std::make_pair(i, j)] = count; ++count;
      }
    }

    Image_neighbor_query neighbor_query(
      image, idx_map, m_image_ptr->get_num_labels(), true);
    Planar_image_region planar_region(
      image, idx_map, m_image_ptr->get_num_labels());
    
    Seed_map seed_map(seeds);
    Region_growing region_growing(
      seeds, neighbor_query, planar_region, seed_map);

    std::vector< std::vector<std::size_t> > regions;
    region_growing.detect(std::back_inserter(regions));

    for (const auto& region : regions) {
      if (region.size() <= 50) {
        const std::size_t new_label = get_best_label(
          neighbor_query, image, region);

        if (new_label == std::size_t(-1))
          continue;
        
        const auto& p = m_image_ptr->get_label_map().at(new_label);
        for (const std::size_t idx : region) {
          const std::size_t ii = image[idx].i;
          const std::size_t jj = image[idx].j;

          auto& cell = original.grid[ii][jj];
          cell.zr = p.x(); cell.zg = p.y(); cell.zb = p.z();
        }
      }
    }

    m_image_ptr->save_image(
      "/Users/monet/Documents/lod/logs/buildings/tmp/image-clean.jpg", original);
  }

  std::size_t get_best_label(
    Image_neighbor_query& neighbor_query,
    const Image& image, 
    const std::vector<std::size_t>& region) {

    const std::size_t ref_label = image[region[0]].label;
    std::vector<std::size_t> nums(m_image_ptr->get_num_labels(), 0);

    if (ref_label == std::size_t(-1))
      return ref_label;

    std::vector<std::size_t> neighbors;
    for (const std::size_t idx : region) {
      neighbors.clear();
      neighbor_query(idx, neighbors);

      for (const std::size_t neighbor : neighbors) {
        if (
          image[neighbor].label != std::size_t(-1) &&
          image[neighbor].label != ref_label) {
          
          nums[image[neighbor].label] += 1;
        }
      }
    }

    std::size_t best_idx = std::size_t(-1);
    std::size_t max_value = 0;
    for (std::size_t i = 0; i < nums.size(); ++i) {
      if (nums[i] >= max_value) {
        max_value = nums[i];
        best_idx = i;
      }
    }
    return best_idx;
  }

  void add_internal_constraints() {

  }

  void create_triangulation() {
    
    auto& tri = m_base.delaunay;
    tri.clear();

    std::vector<Vertex_handle> vhs;
    for (const auto& segment : m_constraints) {
      const auto vh1 = tri.insert(segment.source());
      const auto vh2 = tri.insert(segment.target());

      if (vh1 != vh2)
        tri.insert_constraint(vh1, vh2);
    }

    if (tri.number_of_faces() < 1) {
      m_partition_2.clear(); return;
    }
    compute_visibility(m_base);

    std::map<Face_handle, int> fmap;
    create_faces(m_base, m_partition_2, fmap);
    create_face_neighbors(m_base, fmap, m_partition_2);
    create_edges(m_base, fmap, m_partition_2);
  }

  void compute_visibility(Triangulation& base) const {
    auto& tri = base.delaunay;
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
      } else {
        pface.inside = FT(0); pface.outside = FT(1);
        pface.visibility = Visibility_label::OUTSIDE;
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
};

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_INTERNAL_PARTITION_BUILDER_FROM_IMAGE_2_H
