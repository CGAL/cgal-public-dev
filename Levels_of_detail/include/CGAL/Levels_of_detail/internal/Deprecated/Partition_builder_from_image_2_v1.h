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

#ifndef CGAL_LEVELS_OF_DETAIL_INTERNAL_PARTITION_BUILDER_FROM_IMAGE_2_V1_H
#define CGAL_LEVELS_OF_DETAIL_INTERNAL_PARTITION_BUILDER_FROM_IMAGE_2_V1_H

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
#include <CGAL/Levels_of_detail/internal/Shape_detection/Linear_image_region.h>
#include <CGAL/Levels_of_detail/internal/Spatial_search/Image_neighbor_query.h>
#include <CGAL/Levels_of_detail/internal/Buildings/Building_walls_creator.h>

// Testing.
#include "../../../../../test/Levels_of_detail/include/Saver.h"

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

template<
typename GeomTraits,
typename ImagePointer>
class Partition_builder_from_image_2_v1 {

public:
  using Traits = GeomTraits;

  using FT = typename Traits::FT;
  using Point_2 = typename Traits::Point_2;
  using Point_3 = typename Traits::Point_3;
  using Vector_2 = typename Traits::Vector_2;
  using Vector_3 = typename Traits::Vector_3;
  using Line_2 = typename Traits::Line_2;
  using Line_3 = typename Traits::Line_3;
  using Segment_2 = typename Traits::Segment_2;
  using Plane_3 = typename Traits::Plane_3;
  using Intersect_3 = typename Traits::Intersect_3;

  using Triangulation = internal::Triangulation<Traits>;
  using Face_handle = typename Triangulation::Delaunay::Face_handle;
  using Vertex_handle = typename Triangulation::Delaunay::Vertex_handle;
  using Partition_edge_2 = internal::Partition_edge_2<Traits>;
  using Partition_face_2 = internal::Partition_face_2<Traits>;

  using Partition_2 = internal::Partition_2<Traits>;
  using LF_circulator = typename Triangulation::Delaunay::Line_face_circulator;
  using Size_pair = std::pair<std::size_t, std::size_t>;

  using Building_walls_creator = internal::Building_walls_creator<Traits>;

  struct Pixel {
    Pixel() { }
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
    std::size_t priority = 2;
    std::size_t new_index = std::size_t(-1);
    bool used = false;
  };

  using Image = std::vector<Pixel>;
  using Idx_map = std::map<Size_pair, std::size_t>;
  using Seeds = std::vector<std::size_t>;
  using Seed_map = internal::Seed_property_map;

  using Image_neighbor_query = internal::Image_neighbor_query<Traits, Pixel>;
  using Planar_image_region = internal::Planar_image_region<Traits, Pixel>;
  using Region_growing = internal::Region_growing<
    Seeds, Image_neighbor_query, Planar_image_region, Seed_map>;

  struct Constraint {
    Constraint(
      const Segment_2& segment_,
      const Size_pair& labels_,
      const std::size_t comp_idx_) :
    segment(segment_),
    labels(labels_),
    comp_idx(comp_idx_)
    { }

    const Segment_2 segment;
    const Size_pair labels;
    const std::size_t comp_idx;
  };

  Partition_builder_from_image_2_v1(
    const std::vector<Segment_2>& boundary,
    ImagePointer& image_ptr,
    Partition_2& partition_2) :
  m_boundary(boundary),
  m_image_ptr(image_ptr),
  m_partition_2(partition_2),
  m_pi(static_cast<FT>(CGAL_PI))
  { }

  void build() {
    m_partition_2.clear();

    m_constraints.clear();
    add_boundary_constraints();

    Image image;
    Idx_map idx_map;
    clean_image(image, idx_map);

    add_internal_constraints(image, idx_map);
    create_triangulation();

    std::vector<Segment_2> segments;
    segments.reserve(m_constraints.size());
    for (const auto& constr : m_constraints)
      segments.push_back(constr.segment);

    Saver<Traits> saver;
    saver.save_polylines(segments,
      "/Users/monet/Documents/gf/lod/logs/buildings/tmp/constraints");
  }

private:
  const std::vector<Segment_2>& m_boundary;

  ImagePointer& m_image_ptr;
  Partition_2& m_partition_2;

  std::vector<Constraint> m_constraints;
  Triangulation m_base;

  const FT m_pi;

  void add_boundary_constraints() {
    const Size_pair pair = std::make_pair(std::size_t(-1), std::size_t(-1));
    for (const auto& segment : m_boundary)
      m_constraints.push_back(Constraint(segment, pair, std::size_t(-1)));
  }

  void clean_image(
    Image& image, Idx_map& idx_map) {

    image.clear();
    idx_map.clear();
    Seeds seeds;

    auto& original = m_image_ptr->get_image();

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
      image, idx_map, m_image_ptr->get_num_labels());
    Planar_image_region planar_region(
      image, idx_map, m_image_ptr->get_num_labels());

    Seed_map seed_map(seeds);
    Region_growing region_growing(
      seeds, neighbor_query, planar_region, seed_map);

    std::vector< std::vector<std::size_t> > regions;
    region_growing.detect(std::back_inserter(regions));

    std::cout << "Num label regions: " << regions.size() << std::endl;

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
          image[idx].label = new_label;
        }
      }
    }

    m_image_ptr->save_image(
      "/Users/monet/Documents/gf/lod/logs/buildings/tmp/image-clean.jpg", original);
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

  void add_internal_constraints(
    const Image& image,
    const Idx_map& idx_map) {

    std::vector<Size_pair> pairs;
    find_label_pairs(image, idx_map, pairs);
    for (const auto& pair : pairs)
      add_internal_constraint(pair, image, idx_map);
  }

  void find_label_pairs(
    const Image& image,
    const Idx_map& idx_map,
    std::vector<Size_pair>& pairs) {
    pairs.clear();

    const std::size_t num_labels = m_image_ptr->get_num_labels();
    Image_neighbor_query neighbor_query(
      image, idx_map, num_labels);

    std::set<Size_pair> tmp;
    std::vector<std::size_t> neighbors;
    for (const auto& pixel : image) {
      if (pixel.label == std::size_t(-1))
        continue;

      neighbors.clear();
      neighbor_query(pixel.index, neighbors);

      for (const std::size_t neighbor : neighbors) {
        if (
          image[neighbor].label != std::size_t(-1) &&
          image[neighbor].label != pixel.label) {

          tmp.insert(
            std::make_pair(pixel.label, image[neighbor].label));
        }
      }
    }

    for (const auto& item : tmp) {
      bool found = false;
      for (const auto& pair : pairs) {
        if (
          ( pair.first == item.first && pair.second == item.second ) ||
          ( pair.second == item.first && pair.first == item.second )) {

          found = true;
          break;
        }
      }
      if (!found)
        pairs.push_back(item);
    }

    std::cout << "Num label pairs: " << pairs.size() << std::endl;
  }

  void add_internal_constraint(
    const Size_pair& pair,
    const Image& image,
    const Idx_map& idx_map) {

    std::vector<Pixel> pixels;
    get_pixels(pair, image, idx_map, pixels);

    /*
    auto& original = m_image_ptr->get_image();
    for (const auto& pixel : pixels) {
      auto& cell = original.grid[pixel.i][pixel.j];
      cell.zr = 0;
      cell.zg = 0;
      cell.zb = 0;
    }
    m_image_ptr->save_image(
      "/Users/monet/Documents/gf/lod/logs/buildings/tmp/ridges/image-"
      + std::to_string(pair.first) + "-" + std::to_string(pair.second) +
      ".jpg", original); */

    add_path(pair, image, idx_map, pixels);
  }

  void get_pixels(
    const Size_pair& pair,
    const Image& image,
    const Idx_map& idx_map,
    std::vector<Pixel>& pixels) {

    pixels.clear();
    const std::size_t num_labels = m_image_ptr->get_num_labels();
    Image_neighbor_query neighbor_query(
      image, idx_map, num_labels);

    std::vector<std::size_t> neighbors;
    for (auto& pixel : image) {
      if (pixel.label == pair.first) {

        neighbors.clear();
        neighbor_query(pixel.index, neighbors);

        bool found = false;
        for (const std::size_t neighbor : neighbors) {
          if (image[neighbor].label == pair.second) {
            found = true;
            break;
          }
        }
        if (found)
          pixels.push_back(pixel);
      }
    }
  }

  void add_path(
    const Size_pair& pair,
    const Image& image,
    const Idx_map& idx_map,
    const std::vector<Pixel>& pixels) {

    const std::size_t num_labels = m_image_ptr->get_num_labels();
    Image_neighbor_query neighbor_query(
      image, idx_map, num_labels, true, false);

    Seeds seeds(image.size(), std::size_t(-1));
    for (const auto& pixel : pixels)
      seeds[pixel.index] = pixel.index;
    neighbor_query.make_linear(seeds);

    using Linear_image_region = internal::Linear_image_region<Traits, Pixel>;
    using Region_growing = internal::Region_growing<
      Seeds, Image_neighbor_query, Linear_image_region, Seed_map>;

    Linear_image_region linear_region(
      image, idx_map, m_image_ptr->get_num_labels());

    Seed_map seed_map(seeds);
    Region_growing region_growing(
      seeds, neighbor_query, linear_region, seed_map);

    std::vector< std::vector<std::size_t> > regions;
    region_growing.detect(std::back_inserter(regions));

    // std::cout << "Num components: " << regions.size() << std::endl;
    for (std::size_t i = 0; i < regions.size(); ++i)
      handle_region_v2(i, pair, regions[i], image, idx_map);
  }

  void handle_region_v1(
    const std::size_t region_index,
    const Size_pair& pair,
    const std::vector<std::size_t>& region,
    const Image& image,
    const Idx_map& idx_map) {

    // Find pixels.
    std::vector<Pixel> pixels;
    for (const std::size_t idx : region)
      pixels.push_back(image[idx]);

    const std::size_t num_labels = m_image_ptr->get_num_labels();
    Image_neighbor_query neighbor_query(
      image, idx_map, num_labels, false, false);
    set_priorities(image, neighbor_query, pixels);

    std::sort(pixels.begin(), pixels.end(),
    [](const Pixel& a, const Pixel& b) -> bool {
      return a.priority < b.priority;
    });

    std::size_t count = 0;
    Seeds seeds(image.size(), std::size_t(-1));
    for (auto& pixel : pixels) {
      pixel.new_index = count; ++count;
      seeds[pixel.index] = pixel.new_index;
    }

    // Orient.
    std::vector<Segment_2> segments;
    std::vector<Pixel> path;
    for (auto& pixel : pixels) {
      path.clear();
      if (!pixel.used)
        traverse_path(pixel, neighbor_query, seeds, pixels, path);
      else continue;

      if (path.size() < 2)
        continue;

      std::vector<cv::Point> in, out;
      for (std::size_t k = 0; k < path.size(); ++k) {
        const auto point = cv::Point(path[k].i, path[k].j);
        in.push_back(point);
      }

      cv::approxPolyDP(
        cv::Mat(in), out, 0.0, false);

      segments.clear();
      for (std::size_t k = 0; k < out.size() - 1; ++k) {
        const std::size_t kp = k + 1;
        const auto& q1 = out[k];
        const auto& q2 = out[kp];

        const auto s = m_image_ptr->get_point(q1);
        const auto t = m_image_ptr->get_point(q2);

        segments.push_back(Segment_2(s, t));
      }
    }

    handle_segments_v1(pair, segments);
  }

  void traverse_path(
    Pixel& curr,
    const Image_neighbor_query& neighbor_query,
    const Seeds& seeds,
    std::vector<Pixel>& pixels,
    std::vector<Pixel>& path) {

    curr.used = true;
    path.push_back(curr);

    std::vector<size_t> neighbors;
    neighbor_query(curr.index, neighbors);

    for (const std::size_t idx : neighbors) {
      if (seeds[idx] != std::size_t(-1)) {
        if (!pixels[seeds[idx]].used) {
          auto& next = pixels[seeds[idx]];
          traverse_path(next, neighbor_query, seeds, pixels, path);
        }
      }
    }
  }

  void handle_segments_v1(
    const Size_pair& pair,
    std::vector<Segment_2>& polyline) {

    std::cout << polyline.size() << std::endl;
    for (const auto& segment : polyline)
      m_constraints.push_back(segment);
  }

  void handle_segments_v2(
    const Size_pair& pair,
    std::vector<Segment_2>& polyline) {

    std::cout << "Num segments "
    << pair.first << "-" << pair.second << ": " << polyline.size() << std::endl;

    const auto& plane1 = m_image_ptr->get_plane_map().at(pair.first);
    const auto& plane2 = m_image_ptr->get_plane_map().at(pair.second);

    typename CGAL::cpp11::result_of<
    Intersect_3(Plane_3, Plane_3)>::type result
      = CGAL::intersection(plane1, plane2);

    Line_3 line; bool found = false;
    if (result) {
      if (const Line_3* l = boost::get<Line_3>(&*result)) {
        found = true;
        line = *l;
      }
    }

    if (!found) {
      for (const auto& segment : polyline)
        m_constraints.push_back(segment);
      return;
    }

    const auto v = line.to_vector();
    const auto qq1 = line.point(0);
    const auto qq2 = line.point(1);
    const auto q1 = Point_2(qq1.x(), qq1.y());
    const auto q2 = Point_2(qq2.x(), qq2.y());
    const Line_2 line2 = Line_2(q1, q2);

    auto longest = Segment_2(q1, q2);
    for (auto& segment : polyline) {

      const FT angle = angle_degree_2(longest, segment);
      const FT angle_2 = get_angle_2(angle);

      if (CGAL::abs(angle_2) < FT(10)) {
        const auto p1 = line2.projection(segment.source());
        const auto p2 = line2.projection(segment.target());
        m_constraints.push_back(Segment_2(p1, p2));
      } else {
        m_constraints.push_back(segment);
      }
    }
  }

  FT angle_degree_2(
    const Segment_2& longest, const Segment_2& segment) {

    const Vector_2 v1 =  segment.to_vector();
    const Vector_2 v2 = -longest.to_vector();

    const FT det = CGAL::determinant(v1, v2);
    const FT dot = CGAL::scalar_product(v1, v2);
    const FT angle_rad = static_cast<FT>(
      std::atan2(CGAL::to_double(det), CGAL::to_double(dot)));
    const FT angle_deg = angle_rad * FT(180) / m_pi;
    return angle_deg;
  }

  FT get_angle_2(const FT angle) {

    FT angle_2 = angle;
    if (angle_2 > FT(90)) angle_2 = FT(180) - angle_2;
    else if (angle_2 < -FT(90)) angle_2 = FT(180) + angle_2;
    return angle_2;
  }

  void set_priorities(
    const Image& image,
    const Image_neighbor_query& neighbor_query,
    std::vector<Pixel>& pixels) {

    for (auto& pixel : pixels)
      pixel.priority = 2;

    if (m_image_ptr->get_num_labels() == 1)
      return;

    std::vector<std::size_t> neighbors;
    for (auto& pixel : pixels) {

      neighbors.clear();
      neighbor_query(pixel.index, neighbors);

      std::set<std::size_t> ns;
      ns.insert(pixel.label);

      bool found_boundary = false;
      for (const std::size_t neighbor : neighbors) {
        if (image[neighbor].label == std::size_t(-1)) {
          found_boundary = true; break;
        }
        if (image[neighbor].label != pixel.label)
          ns.insert(image[neighbor].label);
      }

      if (found_boundary) {
        pixel.priority = 0; continue;
      }

      if (ns.size() > 2) {
        pixel.priority = 1; continue;
      }
    }
  }

  void handle_region_v2(
    const std::size_t region_index,
    const Size_pair& pair,
    const std::vector<std::size_t>& region,
    const Image& image,
    const Idx_map& idx_map) {

    std::vector<Pixel> pixels;
    pixels.reserve(region.size());
    for (const std::size_t idx : region)
      pixels.push_back(image[idx]);

    /*
    std::vector<Point_3> points;
    points.reserve(pixels.size());
    for (const auto& pixel : pixels) {
      const auto q = m_image_ptr->get_point(cv::Point(pixel.i, pixel.j));
      points.push_back(Point_3(q.x(), q.y(), FT(0)));
    }

    const CGAL::Color color = CGAL::Color(0, 0, 0);
    const std::string name =
    "/Users/monet/Documents/gf/lod/logs/buildings/tmp/ridges/points-" +
    std::to_string(pair.first) + "-" +
    std::to_string(pair.second) + "-" +
    std::to_string(region_index);

    Saver<Traits> saver;
    saver.export_points(points, color, name); */

    std::vector<Point_2> points;
    points.reserve(pixels.size());
    for (const auto& pixel : pixels) {
      const auto q = m_image_ptr->get_point(cv::Point(pixel.i, pixel.j));
      points.push_back(q);
    }

    std::vector< std::vector<std::size_t> > wall_points_2;
    std::vector<Segment_2> approximate_boundaries_2;

    Building_walls_creator creator(points);
    creator.create_wall_regions(
      0.5,
      0.25,
      25.0,
      0.0,
      wall_points_2);

    creator.create_boundaries(
      wall_points_2,
      approximate_boundaries_2);

    for (const auto& segment : approximate_boundaries_2)
      m_constraints.push_back(
        Constraint(segment, pair, region_index));
  }

  void create_triangulation() {

    auto& tri = m_base.delaunay;
    tri.clear();

    std::vector<Vertex_handle> vhs;
    for (const auto& constr : m_constraints) {
      const auto vh1 = tri.insert(constr.segment.source());
      const auto vh2 = tri.insert(constr.segment.target());

      if (vh1 != vh2) {
        vh1->info().labels.push_back(constr.labels);
        vh2->info().labels.push_back(constr.labels);
        vh1->info().components.push_back(constr.comp_idx);
        vh2->info().components.push_back(constr.comp_idx);
        tri.insert_constraint(vh1, vh2);
      }
    }

    if (tri.number_of_faces() < 1) {
      m_partition_2.clear(); return;
    }
    compute_visibility(m_base);
    clean_triangulation(m_base);

    std::map<Face_handle, int> fmap;
    create_faces(m_base, m_partition_2, fmap);
    create_face_neighbors(m_base, fmap, m_partition_2);
    create_edges(m_base, fmap, m_partition_2);
  }

  void compute_visibility(Triangulation& base) const {

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

            bool found1 = false, found2 = false;
            for (const auto& comp : vh1->info().components) {
              if (comp == std::size_t(-1)) {
                found1 = true; break;
              }
            }
            for (const auto& comp : vh2->info().components) {
              if (comp == std::size_t(-1)) {
                found2 = true; break;
              }
            }

            if (found1 && found2)
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

  void clean_triangulation(Triangulation& base) const {

  }
};

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_INTERNAL_PARTITION_BUILDER_FROM_IMAGE_2_V1_H
