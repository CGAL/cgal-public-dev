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
#include <CGAL/Levels_of_detail/internal/Spatial_search/Image_neighbor_query.h>

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
  using Intersect_3 = typename Traits::Intersect_3;

  using Triangulation    = internal::Triangulation<Traits>;
  using Face_handle      = typename Triangulation::Delaunay::Face_handle;
  using Vertex_handle    = typename Triangulation::Delaunay::Vertex_handle;
  using Partition_edge_2 = internal::Partition_edge_2<Traits>;
  using Partition_face_2 = internal::Partition_face_2<Traits>;

  using Partition_2   = internal::Partition_2<Traits>;
  using LF_circulator = typename Triangulation::Delaunay::Line_face_circulator;
  using Size_pair     = std::pair<std::size_t, std::size_t>;
  using Indices       = std::vector<std::size_t>;
  using Idx_map       = std::map<Size_pair, std::size_t>;
  using Vh_pair       = std::pair<Vertex_handle, Vertex_handle>;

  using Saver = Saver<Traits>;

  struct Pixel {

    std::size_t index = std::size_t(-1);
    std::size_t i = std::size_t(-1), j = std::size_t(-1);
    bool is_boundary = false;
    std::size_t label = std::size_t(-1);
    Point_2 point;
    std::vector<Point_2> duals;
    std::size_t ridge_index = std::size_t(-1);
    Vertex_handle vh;
  };

  struct Image {

    Idx_map idx_map;
    std::vector<Pixel> pixels;
    Indices seeds;

    void clear() {
      idx_map.clear();
      pixels.clear();
      seeds.clear();
    }
  };

  struct Ridge {
    std::vector<Pixel> pixels;
    
    void clear() {
      pixels.clear();
    }
  };

  struct Constraint {
    bool is_boundary = false;
  };

  using Seed_map             = internal::Seed_property_map;
  using Image_neighbor_query = internal::Image_neighbor_query<Traits, Pixel>;
  using Linear_image_region  = internal::Linear_image_region<Traits, Pixel>;
  using Planar_image_region  = internal::Planar_image_region<Traits, Pixel>;

  Partition_builder_from_image_2(
    const std::vector<Segment_2>& boundary,
    ImagePointer& image_ptr,
    Partition_2& partition_2) :
  m_boundary(boundary),
  m_image_ptr(image_ptr),
  m_partition_2(partition_2),
  m_pi(static_cast<FT>(CGAL_PI)) { 

    m_partition_2.clear();
    create_image(m_image);
  }

  void build() {

    clean_image(
      m_image);
    find_label_pairs(
      m_image, m_label_pairs);
    find_ridges(
      m_image, m_label_pairs, m_ridges);
    create_triangulation(
      m_boundary, m_ridges, m_outer_constraints, m_base);

    transform(m_base, m_partition_2);
    save_partition_2(
      "/Users/monet/Documents/lod/logs/buildings/tmp/partition_step_0", false);
    save_ridge_pixels(
      "/Users/monet/Documents/lod/logs/buildings/tmp/ridge_pixels");
  }

  void add_inner_constraints() {

  }

  void compute_visibility() {

    apply_ray_shooting_visibility(
      m_outer_constraints, m_base);
    transform(m_base, m_partition_2);
    save_partition_2(
      "/Users/monet/Documents/lod/logs/buildings/tmp/partition_step_2", false);
  }

  void label_faces() {

    apply_naive_labeling(m_image, m_base);
    transform(m_base, m_partition_2);
    save_partition_2(
      "/Users/monet/Documents/lod/logs/buildings/tmp/partition_step_3", true);
  }

  void optimize() {

  }

private:
  const std::vector<Segment_2>& m_boundary;
  ImagePointer& m_image_ptr;
  Partition_2& m_partition_2;
  const FT m_pi;

  Image m_image;
  Triangulation m_base;
  Saver m_saver;
  
  std::vector<Size_pair> m_label_pairs;
  std::vector<Ridge> m_ridges;
  
  std::map<Vh_pair, Constraint> m_inner_constraints;
  std::map<Vh_pair, Constraint> m_outer_constraints;

  void create_image(Image& image) {
    image.clear();

    std::size_t index = 0;
    auto& original = m_image_ptr->get_image();
    
    Pixel pixel;
    for (std::size_t i = 0; i < original.rows; ++i) {
      for (std::size_t j = 0; j < original.cols; ++j) {

        if (
          i == 0 || j == 0 || 
          i == original.rows - 1 || j == original.cols - 1) {
          
          pixel.label = std::size_t(-1);
          pixel.is_boundary = true;

        } else {
          const auto& cell = original.grid[i][j];
          const std::size_t label = m_image_ptr->get_label(
            cell.zr, cell.zg, cell.zb);
          
          if (label != m_image_ptr->get_num_labels())
            pixel.label = label;
          else 
            pixel.label = std::size_t(-1);
        }
        
        pixel.index = index;
        pixel.i = i; pixel.j = j;
        pixel.point = m_image_ptr->get_point(i, j);

        image.pixels.push_back(pixel);
        
        if (pixel.label == std::size_t(-1))
          image.seeds.push_back(pixel.label);
        else
          image.seeds.push_back(pixel.index);

        image.idx_map[std::make_pair(i, j)] = pixel.index; ++index;
      }
    }
  }

  void clean_image(Image& image) {

    using Region_growing = internal::Region_growing<
      Indices, Image_neighbor_query, Planar_image_region, Seed_map>;

    Image_neighbor_query neighbor_query(
      image.pixels, image.idx_map, false);
    neighbor_query.use_version_4();

    Planar_image_region planar_region(
      image.pixels, image.idx_map);

    Seed_map seed_map(image.seeds);
    Region_growing region_growing(
      image.seeds, neighbor_query, planar_region, seed_map);

    std::vector< std::vector<std::size_t> > regions;
    region_growing.detect(std::back_inserter(regions));

    std::cout << "Num labeled regions: " << regions.size() << std::endl;

    auto& original = m_image_ptr->get_image();
    for (const auto& region : regions) {
      if (region.size() <= 50) {
        const std::size_t new_label = get_best_label(
          image, neighbor_query, region);

        if (new_label == std::size_t(-1))
          continue;
        
        const auto& p = m_image_ptr->get_label_map().at(new_label);
        for (const std::size_t idx : region) {
          const std::size_t i = image.pixels[idx].i;
          const std::size_t j = image.pixels[idx].j;

          auto& cell = original.grid[i][j];
          cell.zr = p.x(); cell.zg = p.y(); cell.zb = p.z();
          image.pixels[idx].label = new_label;
        }
      }
    }

    m_image_ptr->save_image(
      "/Users/monet/Documents/lod/logs/buildings/tmp/image-clean.jpg", original);
  }

  std::size_t get_best_label(
    const Image& image, 
    Image_neighbor_query& neighbor_query,
    const std::vector<std::size_t>& region) {

    const std::size_t ref_label = image.pixels[region[0]].label;
    std::vector<std::size_t> nums(m_image_ptr->get_num_labels(), 0);

    if (ref_label == std::size_t(-1))
      return ref_label;

    std::vector<std::size_t> neighbors;
    for (const std::size_t idx : region) {
      neighbors.clear();
      neighbor_query(idx, neighbors);

      for (const std::size_t neighbor : neighbors) {
        if (
          image.pixels[neighbor].label != std::size_t(-1) &&
          image.pixels[neighbor].label != ref_label) {
          
          nums[image.pixels[neighbor].label] += 1;
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

  void find_label_pairs(
    const Image& image, 
    std::vector<Size_pair>& label_pairs) {
    
    Image_neighbor_query neighbor_query(
      image.pixels, image.idx_map, false);
    neighbor_query.use_version_4();

    std::set<Size_pair> unique;
    std::vector<std::size_t> neighbors;
    for (const auto& pixel : image.pixels) {
      if (pixel.label == std::size_t(-1))
        continue;

      neighbors.clear();
      neighbor_query(pixel.index, neighbors);

      for (const std::size_t neighbor : neighbors) {
        if (
          image.pixels[neighbor].label != std::size_t(-1) &&
          image.pixels[neighbor].label != pixel.label) {

          const std::size_t l1 = pixel.label;
          const std::size_t l2 = image.pixels[neighbor].label;

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
    
    std::cout << "Num label pairs: " << label_pairs.size() << std::endl;
  }

  void find_ridges(
    const Image& image, 
    const std::vector<Size_pair>& label_pairs,
    std::vector<Ridge>& ridges) {

    ridges.clear();
    for (const auto& label_pair : label_pairs)
      extract_ridges(image, label_pair, ridges);
    
    std::cout << "Num ridges: " << ridges.size() << std::endl;
  }

  void extract_ridges(
    const Image& image, 
    const Size_pair& label_pair,
    std::vector<Ridge>& ridges) {

    std::vector<Pixel> rpixels;
    get_ridge_pixels(
      image, label_pair, rpixels);
    add_ridges(
      image, label_pair, rpixels, ridges);
  }

  void get_ridge_pixels(
    const Image& image,
    const Size_pair& label_pair,
    std::vector<Pixel>& rpixels) {

    rpixels.clear();
    Image_neighbor_query neighbor_query(
      image.pixels, image.idx_map, true);
    neighbor_query.use_version_8();

    std::vector<std::size_t> neighbors;
    for (const auto& pixel : image.pixels) {
      if (pixel.label == label_pair.first) {
        const auto& p = pixel.point;
        
        neighbors.clear();
        neighbor_query(pixel.index, neighbors);

        bool found = false;
        for (const std::size_t neighbor : neighbors) {
          const auto& npixel = image.pixels[neighbor];
          const std::size_t nlabel = npixel.label;

          if (nlabel == label_pair.second) {
            found = true; break;
          }
        }

        if (found) {
          
          Pixel rpixel = pixel;
          for (const std::size_t neighbor : neighbors) {
            const auto& npixel = image.pixels[neighbor];
            const std::size_t nlabel = npixel.label;

            if (nlabel != rpixel.label) {
              const auto q = m_image_ptr->get_point(npixel.i, npixel.j);
              const auto d = internal::middle_point_2(p, q);
              rpixel.duals.push_back(d);
            }
          }
          rpixels.push_back(rpixel);
        }
      }
    }
    /* save_ridge(label_pair, 0, rpixels); */
  }

  void add_ridges(
    const Image& image, 
    const Size_pair& label_pair,
    const std::vector<Pixel>& rpixels,
    std::vector<Ridge>& ridges) {

    std::map<std::size_t, std::size_t> mapping;
    Indices seeds(image.pixels.size(), std::size_t(-1));

    for (std::size_t i = 0; i < rpixels.size(); ++i) {
      const std::size_t ridx = rpixels[i].index;
      mapping[ridx] = i;
      seeds[ridx] = ridx;
    }

    Image_neighbor_query neighbor_query(
      image.pixels, image.idx_map, false);
    neighbor_query.use_version_8();
    neighbor_query.use_seeds(seeds);

    using Region_growing = internal::Region_growing<
      Indices, Image_neighbor_query, Linear_image_region, Seed_map>;

    Linear_image_region linear_region(
      image.pixels, image.idx_map);

    Seed_map seed_map(seeds);
    Region_growing region_growing(
      seeds, neighbor_query, linear_region, seed_map);

    std::vector< std::vector<std::size_t> > regions;
    region_growing.detect(std::back_inserter(regions));

    for (std::size_t i = 0; i < regions.size(); ++i) 
      add_ridge(
        label_pair, i, regions[i], mapping, rpixels, ridges);
  }

  void add_ridge(
    const Size_pair& label_pair,
    const std::size_t ridge_index,
    const Indices& region,
    const std::map<std::size_t, std::size_t>& mapping,
    const std::vector<Pixel>& rpixels,
    std::vector<Ridge>& ridges) {

    Ridge ridge; ridge.pixels.reserve(region.size());
    for (const std::size_t idx : region) {
      const std::size_t ridx = mapping.at(idx);
      const auto& rpixel = rpixels[ridx];
      ridge.pixels.push_back(rpixel);
    }
    ridges.push_back(ridge);

    /* save_ridge(label_pair, ridge_index, ridge.pixels); */
  }

  void save_ridge(
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

  void create_triangulation(
    const std::vector<Segment_2>& boundary,
    std::vector<Ridge>& ridges,
    std::map<Vh_pair, Constraint>& outer_constraints,
    Triangulation& base) {
    
    auto& tri = base.delaunay;
    tri.clear();

    Constraint outer_constraint;
    for (const auto& segment : boundary) {
      const auto vh1 = tri.insert(segment.source());
      const auto vh2 = tri.insert(segment.target());

      if (vh1 != vh2) {
        outer_constraint.is_boundary = true;
        tri.insert_constraint(vh1, vh2);
        outer_constraints[std::make_pair(vh1, vh2)] = outer_constraint;
      }
    }

    for (std::size_t i = 0; i < ridges.size(); ++i) {
      for (auto& pixel : ridges[i].pixels) {
        const auto vh = tri.insert(pixel.point);
        vh->info().object_index = pixel.index;
        
        pixel.vh = vh; 
        pixel.ridge_index = i;
      }
    }
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

    m_saver.clear();
    m_saver.export_polygon_soup(vertices, faces, fcolors, path);
  }

  void save_ridge_pixels(
    const std::string path) {

    std::vector<Point_3> pts;
    std::vector< std::vector<Point_3> > points;

    points.reserve(m_ridges.size());
    for (const auto& ridge : m_ridges) {

      pts.clear();
      for (const auto& pixel : ridge.pixels) {
        
        const auto& p = pixel.point;
        const auto& duals = pixel.duals;
        CGAL_assertion(duals.size() != 0);

        Point_2 q;
        internal::compute_barycenter_2(duals, q);
        const auto pt = internal::middle_point_2(p, q);
        pts.push_back(Point_3(pt.x(), pt.y(), FT(0)));
      }
      points.push_back(pts);
    }

    m_saver.clear();
    m_saver.export_points(points, path);
  }

  void apply_ray_shooting_visibility(
    const std::map<Vh_pair, Constraint>& constraints,
    Triangulation& base) {

    auto& tri = base.delaunay;

    // Bbox;
    std::vector<Point_2> points;
    for (auto fit = tri.finite_faces_begin();
    fit != tri.finite_faces_end(); ++fit) {
      const Face_handle fh = static_cast<Face_handle>(fit);
      const auto& p0 = fh->vertex(0)->point();
      const auto& p1 = fh->vertex(1)->point();
      const auto& p2 = fh->vertex(2)->point();
      points.push_back(p0);
      points.push_back(p1);
      points.push_back(p2);
    }

    std::vector<Point_2> bbox;
    CGAL::Identity_property_map<Point_2> pmap;
    internal::bounding_box_2(points, pmap, bbox);
    
    // Visibility.
    for (auto fit = tri.finite_faces_begin();
    fit != tri.finite_faces_end(); ++fit) {
      const Face_handle fh = static_cast<Face_handle>(fit);

      const auto& p0 = fh->vertex(0)->point();
      const auto& p1 = fh->vertex(1)->point();
      const auto& p2 = fh->vertex(2)->point();

      const FT x = (p0.x() + p1.x() + p2.x()) / FT(3);
      const FT y = (p0.y() + p1.y() + p2.y()) / FT(3);
      const Point_2 p = Point_2(x, y);

      FT in = FT(1); FT out = FT(1);
      for (std::size_t i = 0; i < bbox.size(); ++i) {
        const auto& q = bbox[i];
        if (tri.oriented_side(fh, p) == CGAL::ON_NEGATIVE_SIDE)
          continue;

        LF_circulator circ = tri.line_walk(p, q, fh);
        const LF_circulator end = circ;
        if (circ.is_empty()) continue;

        std::size_t inter = 0;
        do {

          LF_circulator f1 = circ; ++circ;
          LF_circulator f2 = circ;

          const bool success = are_neighbors(f1, f2);
          if (!success) break;

          const std::size_t idx = f1->index(f2);
          const auto edge = std::make_pair(f1, idx);
          if (tri.is_constrained(edge)) {
            
            const auto vh1 = f1->vertex( (idx + 1) % 3);
            const auto vh2 = f1->vertex( (idx + 2) % 3);

            const bool found1 = ( 
              constraints.find(std::make_pair(vh1, vh2)) != 
              constraints.end() );
            const bool found2 = ( 
              constraints.find(std::make_pair(vh2, vh1)) != 
              constraints.end() );

            if (found1 || found2)
              ++inter;
          }
          if (tri.is_infinite(f2)) break;
        } while (circ != end);
        if (inter % 2 == 0) out += FT(1);
        else in += FT(1);
      }

      const FT sum = in + out;
      in /= sum; out /= sum;

      if (in > FT(1) / FT(2)) fh->info().interior = true;
      else fh->info().interior = false;
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
        fh->info().label = best_label;
      }
    }
  }
};

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_INTERNAL_PARTITION_BUILDER_FROM_IMAGE_2_H
