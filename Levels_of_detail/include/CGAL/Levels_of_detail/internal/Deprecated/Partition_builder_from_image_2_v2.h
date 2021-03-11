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

#ifndef CGAL_LEVELS_OF_DETAIL_INTERNAL_PARTITION_BUILDER_FROM_IMAGE_2_V2_H
#define CGAL_LEVELS_OF_DETAIL_INTERNAL_PARTITION_BUILDER_FROM_IMAGE_2_V2_H

// STL includes.
#include <map>
#include <set>
#include <vector>
#include <utility>
#include <memory>

// Boost includes.
#define CGAL_DO_NOT_USE_BOYKOV_KOLMOGOROV_MAXFLOW_SOFTWARE
#include <CGAL/internal/Surface_mesh_segmentation/Alpha_expansion_graph_cut.h>

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
#include <CGAL/Levels_of_detail/internal/Buildings/Building_walls_creator.h>
#include <CGAL/Levels_of_detail/internal/Regularization/Segment_merger.h>
#include <CGAL/Levels_of_detail/internal/Regularization/Segment_regularizer.h>

// Testing.
#include "../../../../../test/Levels_of_detail/include/Saver.h"
#include "../../../../../test/Levels_of_detail/include/Utilities.h"

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

template<
typename GeomTraits,
typename ImagePointer>
class Partition_builder_from_image_2_v2 {

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
  using Building_walls_creator = internal::Building_walls_creator<Traits>;
  using Alpha_expansion = CGAL::internal::Alpha_expansion_graph_cut_boost;
  using Segment_merger = internal::Segment_merger<Traits>;
  using Segment_regularizer = internal::Segment_regularizer<Traits>;

  struct Pixel {

    std::size_t index = std::size_t(-1);
    std::size_t i = std::size_t(-1), j = std::size_t(-1);
    bool is_boundary = false;
    std::size_t label = std::size_t(-1);
    Point_2 point;
    std::vector<Point_2> duals;

    std::size_t own_index = std::size_t(-1);
    std::size_t ridge_index = std::size_t(-1);
    Vertex_handle vh;

    Point_2 get_dual_point() const {

      const auto& p  = point;
      const auto& ds = duals;
      CGAL_assertion(ds.size() != 0);

      Point_2 q;
      internal::compute_barycenter_2(ds, q);
      return internal::middle_point_2(p, q);
    }
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
    Size_pair label_pair;

    void clear() {
      pixels.clear();
    }
  };

  struct Constraint {
    bool is_boundary = false;
  };

  struct Item {
    Point_2 point;
  };

  using Seed_map             = internal::Seed_property_map;
  using Image_neighbor_query = internal::Image_neighbor_query<Traits, Pixel>;
  using Linear_image_region  = internal::Linear_image_region<Traits, Pixel>;
  using Planar_image_region  = internal::Planar_image_region<Traits, Pixel>;

  Partition_builder_from_image_2_v2(
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
  m_pi(static_cast<FT>(CGAL_PI)),
  m_simplify(true) {

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
      "/Users/monet/Documents/gf/lod/logs/buildings/tmp/partition_step_0", false);
  }

  void add_constraints() {

    add_inner_constraints(
      m_image, m_ridges, m_inner_constraints, m_base);

    transform(m_base, m_partition_2);
    save_partition_2(
      "/Users/monet/Documents/gf/lod/logs/buildings/tmp/partition_step_1", false);
    save_inner_constraints(
      "/Users/monet/Documents/gf/lod/logs/buildings/tmp/inner_constraints");
    save_ridge_pixels(
      "/Users/monet/Documents/gf/lod/logs/buildings/tmp/ridge_pixels");
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

  void optimize() {

    return;

    apply_gc_optimization(
      m_base, m_optimized_points);
    create_triangulation(
      m_optimized_points, m_base);

    transform(m_base, m_partition_2);
    save_partition_2(
      "/Users/monet/Documents/gf/lod/logs/buildings/tmp/partition_step_4", true);
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
  const bool m_simplify;

  Image m_image;
  Triangulation m_base;
  Saver m_saver;

  std::vector<Size_pair> m_label_pairs;
  std::vector<Ridge> m_ridges;

  std::map<Vh_pair, Constraint> m_inner_constraints;
  std::map<Vh_pair, Constraint> m_outer_constraints;

  std::vector<Item> m_optimized_points;

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
      image.pixels);

    Seed_map seed_map(image.seeds);
    Region_growing region_growing(
      image.seeds, neighbor_query, planar_region, seed_map);

    std::vector<Indices> regions;
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
      "/Users/monet/Documents/gf/lod/logs/buildings/tmp/image-clean.jpg", original);
  }

  std::size_t get_best_label(
    const Image& image,
    Image_neighbor_query& neighbor_query,
    const Indices& region) {

    const std::size_t ref_label = image.pixels[region[0]].label;
    Indices nums(m_image_ptr->get_num_labels(), 0);

    if (ref_label == std::size_t(-1))
      return ref_label;

    Indices neighbors;
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
    Indices neighbors;
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

    Indices neighbors;
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
          FT cx = FT(0), cy = FT(0), ccount = FT(0);
          std::vector<Pixel> ds;
          bool is_corner = false;

          for (const std::size_t neighbor : neighbors) {
            const auto& npixel = image.pixels[neighbor];

            const auto& q = npixel.point;
            const std::size_t nlabel = npixel.label;

            if (nlabel == rpixel.label) // skip this one
              continue;

            if (nlabel == std::size_t(-1)) {
              cx += q.x(); cy += q.y(); ccount += FT(1); // boundary
            } else {
              ds.push_back(npixel); // interior
            }
          }

          // boundary
          if (ccount > FT(0)) {
            cx /= ccount; cy /= ccount;

            const FT lambda = FT(10);
            const FT x = -lambda * p.x() + (FT(1) + lambda) * cx;
            const FT y = -lambda * p.y() + (FT(1) + lambda) * cy;
            rpixel.duals.push_back(Point_2(x, y));
          }

          // interior
          if (ds.size() > 0) {
            if (ds.size() > 1) {

              std::map<std::size_t, bool> nums;
              for (const auto& px : ds)
                nums[px.label] = true;

              if (nums.size() > 1)
                is_corner = true; // corner
            }
            for (const auto& d : ds)
              rpixel.duals.push_back(d.point); // edge
          }

          // add pixel
          rpixels.push_back(rpixel);

          // handle corner
          if (is_corner && m_simplify) {
            for (const auto& d : ds) {
              Pixel corner = d;
              corner.label = rpixel.label;
              corner.duals.clear();
              corner.duals.push_back(rpixel.point);
              rpixels.push_back(corner);
            }
          }
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
    neighbor_query.use_version_4();
    neighbor_query.use_seeds(seeds);

    using Region_growing = internal::Region_growing<
      Indices, Image_neighbor_query, Linear_image_region, Seed_map>;

    Linear_image_region linear_region;
    Seed_map seed_map(seeds);
    Region_growing region_growing(
      seeds, neighbor_query, linear_region, seed_map);

    std::vector<Indices> regions;
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
    ridge.label_pair = label_pair;
    ridges.push_back(ridge);

    /* save_ridge(label_pair, ridge_index, ridge.pixels); */
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
  }

  void add_inner_constraints(
    const Image& image,
    std::vector<Ridge>& ridges,
    std::map<Vh_pair, Constraint>& inner_constraints,
    Triangulation& base) {

    inner_constraints.clear();
    for (std::size_t i = 0; i < ridges.size(); ++i) {
      add_ridge_constraints(
        i, image, ridges[i], inner_constraints, base);
    }
    std::cout << "Constraints are added!" << std::endl;
  }

  void add_ridge_constraints(
    const std::size_t ridge_index,
    const Image& image,
    Ridge& ridge,
    std::map<Vh_pair, Constraint>& inner_constraints,
    Triangulation& base) {

    auto& rpixels = ridge.pixels;
    for (std::size_t i = 0; i < rpixels.size(); ++i) {
      rpixels[i].own_index = i;
      rpixels[i].ridge_index = ridge_index;
    }

    if (m_simplify) {
      add_simplified_constraints(
        ridge_index, ridge, inner_constraints, base);
    } else {
      add_all_ridge_constraints(
        image, ridge, inner_constraints, base);
    }
  }

  void add_simplified_constraints(
    const std::size_t ridge_index,
    Ridge& ridge,
    std::map<Vh_pair, Constraint>& inner_constraints,
    Triangulation& base) {

    const auto& rpixels = ridge.pixels;
    const auto& pair = ridge.label_pair;
    std::vector<Point_2> rpoints, wpoints;

    Line_2 rline; std::size_t count = 0;
    Indices rindices;
    const bool success = intersect_labels(pair, rline);
    if (success) {
      for (const auto& rpixel : rpixels) {
        const auto p = rpixel.get_dual_point();
        const auto q = rline.projection(p);

        const FT distance = internal::distance(p, q);
        if (distance < 0.5) {
          rpoints.push_back(q);
          rindices.push_back(count);
          ++count;
        } else
          wpoints.push_back(p);
      }
    }

    if (rpoints.size() >= 2)
      add_roof_constraints(
        ridge_index, rpoints, rindices, rline, inner_constraints, base);
    if (wpoints.size() >= 2)
      add_wall_constraints(
        ridge_index, wpoints, inner_constraints, base);
  }

  bool intersect_labels(
    const Size_pair& pair, Line_2& line_2) {

    const auto& plane_map = m_image_ptr->get_plane_map();

    const auto& plane1 = plane_map.at(pair.first);
    const auto& plane2 = plane_map.at(pair.second);

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

  void add_roof_constraints(
    const std::size_t ridge_index,
    const std::vector<Point_2>& rpoints,
    const Indices& rindices,
    const Line_2& rline,
    std::map<Vh_pair, Constraint>& inner_constraints,
    Triangulation& base) {

    Point_2 s, t;
    internal::boundary_points_on_line_2(
      rpoints, CGAL::Identity_property_map<Point_2>(),
      rindices, rline, s, t);

    auto& tri = base.delaunay;
    Constraint inner_constraint;
    const auto vh1 = tri.insert(s);
    const auto vh2 = tri.insert(t);

    if (vh1 != vh2) {

      vh1->info().object_index = std::size_t(-1);
      vh1->info().ridge_index = ridge_index;

      vh2->info().object_index = std::size_t(-1);
      vh2->info().ridge_index = ridge_index;

      inner_constraint.is_boundary = false;
      tri.insert_constraint(vh1, vh2);
      inner_constraints[std::make_pair(vh1, vh2)] = inner_constraint;
    }
  }

  void add_wall_constraints(
    const std::size_t ridge_index,
    const std::vector<Point_2>& wpoints,
    std::map<Vh_pair, Constraint>& inner_constraints,
    Triangulation& base) {

    std::vector<Segment_2> wsegments;
    create_approximate_wall_segments(wpoints, wsegments);
    regularize_wall_segments(wsegments);

    auto& tri = base.delaunay;
    Constraint inner_constraint;
    for (const auto& segment : wsegments) {
      const auto vh1 = tri.insert(segment.source());
      const auto vh2 = tri.insert(segment.target());

      if (vh1 != vh2) {

        vh1->info().object_index = std::size_t(-1);
        vh1->info().ridge_index = ridge_index;

        vh2->info().object_index = std::size_t(-1);
        vh2->info().ridge_index = ridge_index;

        inner_constraint.is_boundary = false;
        tri.insert_constraint(vh1, vh2);
        inner_constraints[std::make_pair(vh1, vh2)] = inner_constraint;
      }
    }
  }

  void create_approximate_wall_segments(
    const std::vector<Point_2>& wpoints,
    std::vector<Segment_2>& wsegments) {

    std::vector<Indices> regions;
    Building_walls_creator creator(wpoints);
    creator.create_wall_regions(
      0.5,
      0.5,
      25.0,
      0.25,
      regions);

    creator.create_boundaries(
      regions,
      wsegments);
  }

  void regularize_wall_segments(
    std::vector<Segment_2>& wsegments) {

    std::vector< std::vector<Segment_2> > wcontours(wsegments.size());
    for (std::size_t i = 0; i < wsegments.size(); ++i)
      wcontours[i].push_back(wsegments[i]);

    Segment_regularizer regularizer(
      m_min_length_2, m_angle_bound_2);

    regularizer.compute_multiple_directions(
      m_boundary, wcontours);
    regularizer.regularize_contours(wcontours);

    wsegments.clear();
    for (const auto& wcontour : wcontours)
      for (const auto& wsegment : wcontour)
        wsegments.push_back(wsegment);

    Segment_merger merger(m_ordinate_bound_2);
    merger.merge_segments(wsegments);
    merger.snap_segments(m_boundary, wsegments);
  }

  void add_all_ridge_constraints(
    const Image& image,
    Ridge& ridge,
    std::map<Vh_pair, Constraint>& inner_constraints,
    Triangulation& base) {

    auto& tri = base.delaunay;
    auto& rpixels = ridge.pixels;

    std::map<std::size_t, std::size_t> mapping;
    std::vector<bool> used(image.pixels.size(), false);
    Indices seeds(image.pixels.size(), std::size_t(-1));

    for (std::size_t i = 0; i < rpixels.size(); ++i) {
      const std::size_t ridx = rpixels[i].index;

      mapping[ridx] = i;
      seeds[ridx] = ridx;
    }

    Image_neighbor_query neighbor_query(
      image.pixels, image.idx_map, false);
    neighbor_query.use_version_4();
    neighbor_query.use_seeds(seeds);

    Indices neighbors;
    Constraint inner_constraint;

    for (auto& rpixel : rpixels) {
      if (used[rpixel.index])
        continue;

      used[rpixel.index] = true;
      neighbors.clear();
      neighbor_query(rpixel.index, neighbors);

      for (const std::size_t neighbor : neighbors) {
        auto& npixel = rpixels[mapping.at(neighbor)];

        if (
          seeds[npixel.index] != std::size_t(-1) &&
          used[npixel.index] != true) {

          const auto vh1 = tri.insert(rpixel.get_dual_point());
          const auto vh2 = tri.insert(npixel.get_dual_point());

          if (vh1 != vh2) {

            vh1->info().object_index = rpixel.own_index;
            vh1->info().ridge_index = rpixel.ridge_index;
            rpixel.vh = vh1;

            vh2->info().object_index = npixel.own_index;
            vh2->info().ridge_index = npixel.ridge_index;
            npixel.vh = vh2;

            inner_constraint.is_boundary = false;
            tri.insert_constraint(vh1, vh2);
            inner_constraints[std::make_pair(vh1, vh2)] = inner_constraint;
          }
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

  void apply_gc_optimization(
    const Triangulation& base,
    std::vector<Item>& items) {

    items.clear();

    mark_vertices();

    mark_edges();

    create_labels();

    /*
    const FT beta = FT(1);

    std::vector<std::size_t> labels;
    set_initial_labels(base, labels);

    std::vector<Size_pair> edges;
    std::vector<double> edge_weights;
    set_graphcut_edges(beta, base, edges, edge_weights);

    std::vector< std::vector<double> > cost_matrix;
    set_cost_matrix(base, cost_matrix);

    Alpha_expansion graphcut;
    graphcut(edges, edge_weights, cost_matrix, labels); */
  }

  void mark_vertices() {

  }

  void mark_edges() {

  }

  void create_labels() {

  }

  void create_triangulation(
    std::vector<Item>& items,
    Triangulation& base) {

    auto& tri = base.delaunay;
    tri.clear();

    for (const auto& item : items) {
      const auto& point = item.point;
      tri.insert(point);
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
        const auto pt = pixel.get_dual_point();
        pts.push_back(Point_3(pt.x(), pt.y(), FT(0)));
      }
      points.push_back(pts);
    }

    m_saver.clear();
    m_saver.export_points(points, path);
  }

  void save_inner_constraints(
    const std::string path) {

    std::vector<Segment_2> segments;
    segments.reserve(m_inner_constraints.size());

    for (const auto& pair : m_inner_constraints) {
      const auto vh1 = pair.first.first;
      const auto vh2 = pair.first.second;

      const auto& s = vh1->point();
      const auto& t = vh2->point();
      segments.push_back(Segment_2(s, t));
    }

    m_saver.clear();
    m_saver.save_polylines(segments, path);
  }
};

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_INTERNAL_PARTITION_BUILDER_FROM_IMAGE_2_V2_H
