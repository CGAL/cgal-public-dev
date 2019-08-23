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

#ifndef CGAL_LEVELS_OF_DETAIL_INTERNAL_GENERIC_SIMPLIFIER_H
#define CGAL_LEVELS_OF_DETAIL_INTERNAL_GENERIC_SIMPLIFIER_H

// STL includes.
#include <memory>
#include <vector>
#include <utility>
#include <stdio.h>

// CGAL includes.
#include <CGAL/assertions.h>
#include <CGAL/barycenter.h>
#include <CGAL/property_map.h>

#include <CGAL/Random.h>
#include <CGAL/IO/Color.h>

#define CGAL_DO_NOT_USE_BOYKOV_KOLMOGOROV_MAXFLOW_SOFTWARE
#include <CGAL/internal/Surface_mesh_segmentation/Alpha_expansion_graph_cut.h>

// LOD includes.
#include <CGAL/Levels_of_detail/enum.h>

// Internal includes.
#include <CGAL/Levels_of_detail/internal/utils.h>
#include <CGAL/Levels_of_detail/internal/struct.h>

// Testing.
#include "../../../../../test/Levels_of_detail/include/Saver.h"

// Simplification.
#include <CGAL/Levels_of_detail/internal/Simplification/Alpha_shapes_filtering_2.h>

// Spatial search.
#include <CGAL/Levels_of_detail/internal/Spatial_search/K_neighbor_query.h>

// OpenCV.
#include "opencv2/opencv.hpp"
#include "opencv2/core/core.hpp"
#include "opencv2/core/utility.hpp"

#include "opencv2/imgcodecs.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

  template<
  typename GeomTraits,
  typename PointMap3>
  class Generic_simplifier {

  public:
    using Traits = GeomTraits;
    using Point_map_3 = PointMap3;

    using FT = typename Traits::FT;
    using Point_2 = typename Traits::Point_2;
    using Point_3 = typename Traits::Point_3;
    using Vector_2 = typename Traits::Vector_2;
    using Vector_3 = typename Traits::Vector_3;
    using Plane_3 = typename Traits::Plane_3;

    using Points_2 = std::vector<Point_2>;
    using Points_3 = std::vector<Point_3>;
    using Indices = std::vector<std::size_t>;

    using Saver = Saver<Traits>;
    using Color = CGAL::Color;

    using Size_pair = std::pair<std::size_t, std::size_t>;

    using Identity_map_2 = 
    CGAL::Identity_property_map<Point_2>;

    using Alpha_shapes_filtering_2 = 
    internal::Alpha_shapes_filtering_2<Traits>;

    using Alpha_expansion = CGAL::internal::Alpha_expansion_graph_cut_boost;
    using Triangulation = internal::Triangulation<Traits>;
    using Location_type = typename Triangulation::Delaunay::Locate_type;

    using Pair = std::pair<Point_2, FT>;
    using Pair_map = CGAL::First_of_pair_property_map<Pair>;
    using K_neighbor_query =
    internal::K_neighbor_query<Traits, std::vector<Pair>, Pair_map>;

    struct Cluster_item {
      Cluster_item(
        const Point_3 _point, 
        const std::size_t _roof_idx) :
      input_point(_point),
      roof_idx(_roof_idx) { }
      
      Point_3 input_point;
      Point_3 final_point;
      std::size_t roof_idx;
    };

    using Cell_id = std::pair<long, long>;
    using Cell_data = std::vector<std::size_t>;
    using Grid = std::map<Cell_id, Cell_data>;

    struct Image_cell {
      
      Image_cell() : 
      roof_idx(std::size_t(-1)),
      zr(FT(255)), zg(FT(255)), zb(FT(255)),
      is_interior(false) { }

      std::size_t roof_idx;
      FT zr, zg, zb;
      bool is_interior;
    };

    struct Image {

      Image() : rows(0), cols(0) { }
      Image(const std::size_t _rows, const std::size_t _cols) : 
      rows(_rows), cols(_cols) { 
        resize(_rows, _cols);
      }

      void clear() {
        rows = 0;
        cols = 0;
        grid.clear();
      }

      void resize(const std::size_t _rows, const std::size_t _cols) {
        rows = _rows; cols = _cols;
        grid.resize(_rows);
        for (auto& pixels : grid)
          pixels.resize(_cols);
      }

      void create_pixel(
        const std::size_t i, const std::size_t j,
        const std::size_t roof_idx,
        const bool is_interior,
        const FT zr, const FT zg, const FT zb) {

        auto& pixel = grid[i][j];

        pixel.roof_idx = roof_idx;
        pixel.zr = zr;
        pixel.zg = zg;
        pixel.zb = zb;
        pixel.is_interior = is_interior;
      }

      std::size_t rows, cols;
      std::vector< std::vector<Image_cell> > grid;
    };

    struct Pixel {

      Pixel(
        const Point_2& p,
        const std::size_t _i, const std::size_t _j,
        const bool _is_interior) :
      point(Point_3(p.x(), p.y(), FT(0))),
      i(_i), j(_j),
      is_interior(_is_interior) 
      { }

      Point_3 point;
      std::size_t i;
      std::size_t j;
      bool is_interior;
    };

    struct Mycolor {

      Mycolor() : 
      zr(FT(255)), zg(FT(255)), zb(FT(255)) { }
      Mycolor(const FT zr_, const FT zg_, const FT zb_) :
      zr(zr_), zg(zg_), zb(zb_) { }
      FT zr, zg, zb;
    };

    using OpenCVImage = cv::Mat;

    Generic_simplifier(
      const Indices& input_range,
      const Point_map_3 point_map_3,
      const FT grid_cell_width_2,
      const FT alpha_shape_size_2,
      const FT graph_cut_beta_2) :
    m_input_range(input_range),
    m_point_map_3(point_map_3),
    m_grid_cell_width_2(grid_cell_width_2),
    m_alpha_shape_size_2(alpha_shape_size_2),
    m_val_min(+internal::max_value<FT>()),
    m_val_max(-internal::max_value<FT>()),
    m_num_labels(0),
    m_rows_min(+internal::max_value<long>()),
    m_rows_max(-internal::max_value<long>()),
    m_cols_min(+internal::max_value<long>()),
    m_cols_max(-internal::max_value<long>()),
    m_pixels_per_cell(9),
    m_samples_per_face(20),
    m_beta(graph_cut_beta_2),
    m_k(FT(6)) 
    { }

    void create_cluster() {

      m_cluster.clear();
      m_cluster.reserve(m_input_range.size());

      for (const std::size_t idx : m_input_range) {
        const Point_3& point = get(m_point_map_3, idx);
        m_cluster.push_back(Cluster_item(point, 0));
        m_val_min = CGAL::min(point.z(), m_val_min);
        m_val_max = CGAL::max(point.z(), m_val_max);
      }
      m_num_labels = 1;
      save_cluster("/Users/monet/Documents/lod/logs/buildings/tmp/cluster");
    }

    void create_cluster_from_regions(
      const std::vector<Indices>& regions) {

      std::vector<Points_3> roofs;
      create_sampled_roofs(regions, roofs);

      std::size_t num_points = 0;
      for (const auto& roof : roofs)
        num_points += roof.size();

      m_cluster.clear();
      m_cluster.reserve(num_points);

      for (std::size_t i = 0; i < roofs.size(); ++i) {
        for (const auto& point : roofs[i]) {

          m_cluster.push_back(Cluster_item(point, i));
          m_val_min = CGAL::min(point.z(), m_val_min);
          m_val_max = CGAL::max(point.z(), m_val_max);
        }
      }
      m_num_labels = roofs.size();
      save_cluster("/Users/monet/Documents/lod/logs/buildings/tmp/cluster");
    }

    void transform_cluster() {

      std::vector<Point_2> points;
      points.reserve(m_cluster.size());

      for (const auto& item : m_cluster)
        points.push_back(internal::point_2_from_point_3(item.input_point));

      Vector_2 dir;
      internal::estimate_direction_2(points, dir);
      const Vector_2 y_dir = Vector_2(FT(0), FT(1));

			internal::compute_angle_2(dir, y_dir, m_angle_2d);
      internal::compute_barycenter_2(points, m_b);
      
      for (Point_2& p : points)
        internal::rotate_point_2(m_angle_2d, m_b, p);

      CGAL::Identity_property_map<Point_2> pmap;
      std::vector<Point_2> bbox;
      internal::bounding_box_2(points, pmap, bbox);

      m_tr = bbox[0];
      for (Point_2& p : points)
        internal::translate_point_2(m_tr, p);

      for (std::size_t i = 0; i < points.size(); ++i) {
        const Point_2& p = points[i];
        m_cluster[i].final_point = 
          Point_3(p.x(), p.y(), m_cluster[i].input_point.z());
      }
    }

    void create_grid() {
      CGAL_assertion(m_cluster.size() >= 3);

      Cell_id cell_id;
      m_grid.clear();

      for (std::size_t i = 0; i < m_cluster.size(); ++i) {
        const auto& item = m_cluster[i];
        
        const Point_3& point = item.final_point;
        get_cell_id(point, cell_id);
        m_grid[cell_id].push_back(i);
      }
      save_grid("/Users/monet/Documents/lod/logs/buildings/tmp/grid");
    }

    void create_image() {

      const std::size_t rowsdiff = std::size_t(m_rows_max - m_rows_min);
      const std::size_t colsdiff = std::size_t(m_cols_max - m_cols_min);
      const std::size_t rows = (rowsdiff + 3); // +2 = +1 (diff pixel) +2 (boundary pixels)
      const std::size_t cols = (colsdiff + 3); // +2 = +1 (diff pixel) +2 (boundary pixels)

      std::cout << "Resolution (original): " << cols << "x" << rows << std::endl;
      std::cout << "Cols: " << colsdiff << " Rows: " << rowsdiff << std::endl;
      std::cout << "Val min: " << m_val_min << " Val max: " << m_val_max << std::endl;

      m_image.clear();
      m_image.resize(rows, cols);

      initialize_image(m_image);
      save_image("/Users/monet/Documents/lod/logs/buildings/tmp/image-origin.jpg", m_image);
      create_label_map(m_image);

      inpaint_image_opencv(m_image);
      update_interior_pixels(m_image);
      save_point_cloud("/Users/monet/Documents/lod/logs/buildings/tmp/point-cloud", m_image);
      save_image("/Users/monet/Documents/lod/logs/buildings/tmp/image-paints.jpg", m_image);
      
      apply_graphcut(m_image);
      save_image("/Users/monet/Documents/lod/logs/buildings/tmp/image-gcuted.jpg", m_image);
    }

    void get_outer_boundary_points_2(
      Points_2& boundary_points_2) {

      boundary_points_2.clear();

      const Point_2 tr = Point_2(-m_tr.x(), -m_tr.y());
      std::vector<std::size_t> ni, nj;
      for (long i = 1; i < m_image.rows - 1; ++i) {
        for (long j = 1; j < m_image.cols - 1; ++j) {
          get_grid_neighbors_4(i, j, ni, nj);

          for (std::size_t k = 0; k < 4; ++k) {
            const long ii = ni[k];
            const long jj = nj[k];

            if (is_outer_boundary_pixel(i, j, ii, jj)) {

              Point_2 p = get_point_from_id(i, j);
              Point_2 q = get_point_from_id(ii, jj);
              
              internal::translate_point_2(tr, p);
              internal::translate_point_2(tr, q);

              internal::rotate_point_2(-m_angle_2d, m_b, p);
              internal::rotate_point_2(-m_angle_2d, m_b, q);

              boundary_points_2.push_back(internal::middle_point_2(p, q));
            }
          }
        }
      }
    }

    void get_regular_points(
      std::vector< std::pair<Point_2, bool> >& points) {
      
      std::vector<Pixel> point_cloud;
      create_point_cloud(m_image, point_cloud);
      
      points.clear();
      points.reserve(point_cloud.size());
      const Point_2 tr = Point_2(-m_tr.x(), -m_tr.y());

      for (const auto& pixel : point_cloud) {
        Point_2 p = Point_2(pixel.point.x(), pixel.point.y());

        internal::translate_point_2(tr, p);
        internal::rotate_point_2(-m_angle_2d, m_b, p);

        points.push_back(std::make_pair(p, pixel.is_interior));
      }

      save_regular_points(
        points, "/Users/monet/Documents/lod/logs/buildings/tmp/visibility_points");
    }

    void get_interior_points(
      const Triangulation& tri,
      const Indices& cluster,
      std::vector<Point_3>& points) {

      std::vector<Pair> pairs;
      pairs.reserve(cluster.size());
      for (const std::size_t idx : cluster) {
        const auto& p = get(m_point_map_3, idx);
        pairs.push_back(std::make_pair(Point_2(p.x(), p.y()), p.z()));
      }

      Pair_map pmap;
      K_neighbor_query neighbor_query(pairs, m_k, pmap);

      std::vector<Pixel> point_cloud;
      create_point_cloud(m_image, point_cloud);

      points.clear();
      const Point_2 tr = Point_2(-m_tr.x(), -m_tr.y());

      for (const auto& pixel : point_cloud) {
        if (!pixel.is_interior) continue;

        Point_2 p = Point_2(pixel.point.x(), pixel.point.y());

        internal::translate_point_2(tr, p);
        internal::rotate_point_2(-m_angle_2d, m_b, p);

        Location_type type; int stub;
        const auto fh = tri.delaunay.locate(p, type, stub);
        if (
          type == Triangulation::Delaunay::FACE &&
          !tri.delaunay.is_infinite(fh) &&
          fh->info().tagged) {

          const FT height = get_height(p, pairs, neighbor_query);
          points.push_back(Point_3(p.x(), p.y(), height));
        }
      }

      m_saver.export_points(
        points, 
        Color(0, 0, 0), 
        "/Users/monet/Documents/lod/logs/buildings/tmp/better_cluster");
    }

  private:
    const Indices& m_input_range;
    const Point_map_3 m_point_map_3;

    const FT m_grid_cell_width_2;
    const FT m_alpha_shape_size_2;

    // Cluster.
    std::vector<Cluster_item> m_cluster;
    FT m_val_min, m_val_max;
    std::size_t m_num_labels;
    
    // Transform.
    Point_2 m_b, m_tr;
    FT m_angle_2d;
    
    // Grid.
    Grid m_grid;
    long m_rows_min, m_rows_max;
    long m_cols_min, m_cols_max;
    
    // Image.
    Image m_image;
    std::map<std::size_t, Mycolor> m_label_map;
    const std::size_t m_pixels_per_cell;

    const std::size_t m_samples_per_face;
    const FT m_beta;
    const FT m_k;

    Saver m_saver;

    void create_sampled_roofs(
      const std::vector<Indices>& regions,
      std::vector<Points_3>& roofs) {

      roofs.clear();
      roofs.reserve(regions.size());

      Points_3 roof; Plane_3 plane;
      for (const auto& region : regions) {
        roof.clear();
  
        internal::plane_from_points_3(
        m_input_range, m_point_map_3, region, plane);
        internal::project_on_plane_3(
        m_input_range, m_point_map_3, region, plane, roof);
        sample_roof_region(plane, roof);
        roofs.push_back(roof);
      }
    }

    void sample_roof_region(
      const Plane_3& plane, Points_3& roof) {

      Point_3 b;
      internal::compute_barycenter_3(roof, b);

      Points_2 points;
      points.reserve(roof.size());
      for (const auto& p : roof) {
        const Point_2 q = internal::to_2d(p, b, plane);
        points.push_back(q);
      }
      apply_filtering(points);

      roof.clear();
      for (const auto& p : points) {
        const Point_3 q = internal::to_3d(p, b, plane);
        roof.push_back(q);
      }
    }

    void apply_filtering(std::vector<Point_2>& points) {

      const std::size_t nump = points.size();
      Alpha_shapes_filtering_2 filtering(m_alpha_shape_size_2);
      const FT sampling_2 = m_alpha_shape_size_2 / FT(2);

      Identity_map_2 identity_map_2;
      filtering.add_points(points, identity_map_2);
      points.clear(); 
      filtering.get_samples(sampling_2, m_samples_per_face, points);
    }

    FT get_height(
      const Point_2& p,
      const std::vector<Pair>& pairs,
      K_neighbor_query& neighbor_query) {
      
      Indices neighbors;
      neighbor_query(p, neighbors);

      FT avg_height = FT(0);
      for (const std::size_t idx : neighbors)
        avg_height += pairs[idx].second;
      avg_height /= static_cast<FT>(neighbors.size());

      return avg_height;
    }

    void get_cell_id(
      const Point_3& point, 
      Cell_id& cell_id) {

      const long id_x = get_id_value(point.x());
      const long id_y = get_id_value(point.y());
      
      cell_id = std::make_pair(id_x, id_y);

      m_cols_min = CGAL::min(id_x, m_cols_min);
      m_rows_min = CGAL::min(id_y, m_rows_min);

      m_cols_max = CGAL::max(id_x, m_cols_max);
      m_rows_max = CGAL::max(id_y, m_rows_max);
    }

    long get_id_value(const FT value) {

      CGAL_precondition(m_grid_cell_width_2 > FT(0));
      const long id = static_cast<long>(
        CGAL::to_double(value / m_grid_cell_width_2));
      if (value >= FT(0)) return id;
      return id - 1;
    }

    void initialize_image(
      Image& image) {

      std::size_t numcells = 0;
      for (std::size_t i = 1; i < image.rows - 1; ++i) {
        for (std::size_t j = 1; j < image.cols - 1; ++j) {

          const long id_x = get_id_x(j-1);
          const long id_y = get_id_y(i-1);
          
          const Cell_id cell_id = std::make_pair(id_x, id_y);
          if (m_grid.find(cell_id) != m_grid.end()) {
            ++numcells;

            const auto& indices = m_grid.at(cell_id);
            initialize_pixel(i, j, indices, image);
          }
        }
      }
      std::cout << "Num cells: " << m_grid.size() << " : " << numcells << std::endl;
    }

    long get_id_x(const std::size_t j) {
      return m_cols_min + long(j);
    }

    long get_id_y(const std::size_t i) {
      return m_rows_max - long(i);
    }

    void initialize_pixel(
      const std::size_t i, const std::size_t j,
      const Cell_data& indices,
      Image& image) {

      std::size_t roof_idx = std::size_t(-1);
      FT zr = FT(255), zg = FT(255), zb = FT(255);
      get_pixel_data(indices, roof_idx, zr, zg, zb);
      image.create_pixel(i, j, roof_idx, true, zr, zg, zb);
    }

    void get_pixel_data(
      const Cell_data& indices,
      std::size_t& roof_idx,
      FT& zr, FT& zg, FT& zb) {

      std::vector<int> tmp(m_num_labels, 0);
      for (const std::size_t idx : indices)
        tmp[m_cluster[idx].roof_idx] += 1;

      std::size_t final_idx = std::size_t(-1); int max_num = -1;
      for (std::size_t i = 0; i < tmp.size(); ++i) {
        if (tmp[i] > max_num) {
          final_idx = i;
          max_num = tmp[i];
        }
      }

      roof_idx = final_idx;
      Random rand(roof_idx);
      zr = static_cast<FT>(64 + rand.get_int(0, 192));
      zg = static_cast<FT>(64 + rand.get_int(0, 192));
      zb = static_cast<FT>(64 + rand.get_int(0, 192));
    }

    void create_label_map(
      const Image& image) {

      m_label_map.clear();
      for (std::size_t i = 1; i < image.rows - 1; ++i) {
        for (std::size_t j = 1; j < image.cols - 1; ++j) {
          const auto& cell = image.grid[i][j];
          if (cell.roof_idx != std::size_t(-1)) {

            m_label_map[cell.roof_idx] = Mycolor(
              cell.zr, cell.zg, cell.zb);
          }
        }
      }
    }

    void inpaint_image_opencv(
      Image& image) {
      
      OpenCVImage input(
        image.rows, 
        image.cols, 
        CV_8UC3, cv::Scalar(255, 255, 255));

      OpenCVImage mask(
        image.rows, 
        image.cols, 
        CV_8U, cv::Scalar(0, 0, 0));

      for (std::size_t i = 0; i < image.rows; ++i) {
        for (std::size_t j = 0; j < image.cols; ++j) {
          const uchar zr = saturate_z(image.grid[i][j].zr);
          const uchar zg = saturate_z(image.grid[i][j].zg);
          const uchar zb = saturate_z(image.grid[i][j].zb);
          
          cv::Vec3b& bgr = input.at<cv::Vec3b>(i, j);
          bgr[0] = zb;
          bgr[1] = zg;
          bgr[2] = zr;

          if (!image.grid[i][j].is_interior) {
            unsigned char& val = mask.at<unsigned char>(i, j);
            val = static_cast<unsigned char>(255);
          }
        }
      }

      OpenCVImage inpainted;
      inpaint(input, mask, inpainted, 0, cv::INPAINT_TELEA);

      Image colored(image.rows, image.cols);
      for (std::size_t i = 1; i < colored.rows - 1; ++i) {
        for (std::size_t j = 1; j < colored.cols - 1; ++j) {
          const cv::Vec3b& bgr = inpainted.at<cv::Vec3b>(i, j);
          
          const std::size_t roof_idx = image.grid[i][j].roof_idx;
          const bool is_interior = image.grid[i][j].is_interior;

          const FT zr = FT(static_cast<std::size_t>(bgr[2]));
          const FT zg = FT(static_cast<std::size_t>(bgr[1]));
          const FT zb = FT(static_cast<std::size_t>(bgr[0]));

          colored.create_pixel(i, j, roof_idx, is_interior, zr, zg, zb);
        }
      }
      image = colored;
    }

    void update_interior_pixels(Image& image) {

      std::vector<Pixel> point_cloud;
      create_point_cloud(image, point_cloud);

      Points_3 points;
      create_input_points(points);

      CGAL::Identity_property_map<Point_3> pmap;
      Alpha_shapes_filtering_2 filtering(m_alpha_shape_size_2);
      filtering.add_points(points, pmap);
      filtering.set_interior_labels(point_cloud);

      for (const auto& pixel : point_cloud) {
        if (pixel.is_interior) {
          image.grid[pixel.i][pixel.j].is_interior = true;
        } else {
          image.grid[pixel.i][pixel.j].is_interior = false;
          image.grid[pixel.i][pixel.j].zr = FT(255);
          image.grid[pixel.i][pixel.j].zg = FT(255);
          image.grid[pixel.i][pixel.j].zb = FT(255);
          image.grid[pixel.i][pixel.j].roof_idx = std::size_t(-1);
        }
      }
    }

    void create_input_points(Points_3& points) {
      
      points.clear();
      points.reserve(m_input_range.size());
      for (const std::size_t idx : m_input_range) {
        const Point_3& p = get(m_point_map_3, idx);

        Point_2 q = Point_2(p.x(), p.y());
        internal::rotate_point_2(m_angle_2d, m_b, q);
        internal::translate_point_2(m_tr, q);
        points.push_back(Point_3(q.x(), q.y(), FT(0)));
      }
    }

    void create_point_cloud(
      const Image& image,
      std::vector<Pixel>& point_cloud) {

      point_cloud.clear();
      for (std::size_t i = 0; i < image.rows; ++i) {
        for (std::size_t j = 0; j < image.cols; ++j) {
          const auto& cell = image.grid[i][j];
          
          const bool is_interior = cell.is_interior;
          const Point_2 p = get_point_from_id(i, j);
          point_cloud.push_back(Pixel(p, i, j, is_interior));
        }
      }
    }

    Point_2 get_point_from_id(
      const std::size_t i, const std::size_t j) {

      const long id_x = get_id_x(j);
      const long id_y = get_id_y(i);

      const FT x = get_coordinate(id_x);
      const FT y = get_coordinate(id_y);

      return Point_2(x, y);
    }

    const FT get_coordinate(long id) {

      CGAL_precondition(m_grid_cell_width_2 > FT(0));
      if (id < 0) id = id + 1;
      const FT half = m_grid_cell_width_2 / FT(2);
      const FT value = static_cast<FT>(id) * m_grid_cell_width_2 + half;
      return value;
    }

    void apply_graphcut(
      Image& image) {

      std::map<Size_pair, std::size_t> idx_map;
      set_idx_map(image, idx_map);

      std::vector<std::size_t> labels;
      set_initial_labels(image, idx_map, labels);
      apply_new_labels(idx_map, labels, image);

      save_image("/Users/monet/Documents/lod/logs/buildings/tmp/image-labels.jpg", image);

      std::vector<Size_pair> edges;
      std::vector<double> edge_weights;
      set_graphcut_edges(image, idx_map, edges, edge_weights);
      
      std::vector< std::vector<double> > cost_matrix;
      set_cost_matrix(image, idx_map, cost_matrix);

      compute_graphcut(edges, edge_weights, cost_matrix, labels);
      apply_new_labels(idx_map, labels, image);
    }

    void set_idx_map(
      const Image& image,
      std::map<Size_pair, std::size_t>& idx_map) {

      idx_map.clear();
      std::size_t pixel_idx = 0;
      for (std::size_t i = 1; i < image.rows - 1; ++i) {
        for (std::size_t j = 1; j < image.cols - 1; ++j) {
          // if (!image.grid[i][j].is_interior) continue;

          idx_map[std::make_pair(i, j)] = pixel_idx;
          ++pixel_idx;
        }
      }
    }

    void set_initial_labels(
      const Image& image,
      const std::map<Size_pair, std::size_t>& idx_map,
      std::vector<std::size_t>& labels) {

      labels.clear();
      labels.resize(idx_map.size());
      for (std::size_t i = 1; i < image.rows - 1; ++i) {
        for (std::size_t j = 1; j < image.cols - 1; ++j) {
          // if (!image.grid[i][j].is_interior) continue;

          const std::size_t label = get_label(
            image.grid[i][j].zr,
            image.grid[i][j].zg,
            image.grid[i][j].zb);
          const std::size_t pixel_idx = idx_map.at(std::make_pair(i, j));
          labels[pixel_idx] = label;
        }
      }
    }

    std::size_t get_label(
      const FT zr, const FT zg, const FT zb) {

      if (zr == FT(255) && zg == FT(255) && zb == FT(255)) // tmp
        return m_num_labels;

      FT d_max = FT(-1); std::size_t label = std::size_t(-1);
      for (const auto& pair: m_label_map) {
        const FT zr_diff = zr - pair.second.zr;
        const FT zg_diff = zg - pair.second.zg;
        const FT zb_diff = zb - pair.second.zb;

        const double r = CGAL::to_double(zr_diff * zr_diff);
        const double g = CGAL::to_double(zg_diff * zg_diff);
        const double b = CGAL::to_double(zb_diff * zb_diff);
        
        const FT d = static_cast<FT>(CGAL::sqrt(r + g + b));
        if (d > d_max) {
          d_max = d; label = pair.first;
        }
      }
      return label;
    }

    void apply_new_labels(
      const std::map<Size_pair, std::size_t>& idx_map,
      const std::vector<std::size_t>& labels,
      Image& image) {

      Image labeled(image.rows, image.cols);
      for (std::size_t i = 1; i < image.rows - 1; ++i) {
        for (std::size_t j = 1; j < image.cols - 1; ++j) {
          // if (!image.grid[i][j].is_interior) continue;
          
          const std::size_t pixel_idx = idx_map.at(std::make_pair(i, j));
          
          Mycolor color; 
          bool is_interior = image.grid[i][j].is_interior;

          if (labels[pixel_idx] == m_num_labels) {
            color = Mycolor(FT(255), FT(255), FT(255));
            is_interior = false;
          } else {
            color = m_label_map.at(labels[pixel_idx]);
            is_interior = true;
          }
          
          const std::size_t roof_idx = image.grid[i][j].roof_idx;
          labeled.create_pixel(i, j, roof_idx, is_interior, 
            color.zr, color.zg, color.zb);
        }
      }
      image = labeled;
    }

    void set_graphcut_edges(
      const Image& image,
      const std::map<Size_pair, std::size_t>& idx_map,
      std::vector<Size_pair>& edges,
      std::vector<double>& edge_weights) {

      edges.clear();
      edge_weights.clear();
      std::vector<std::size_t> ni, nj;
      for (std::size_t i = 1; i < image.rows - 1; ++i) {
        for (std::size_t j = 1; j < image.cols - 1; ++j) {
          // if (!image.grid[i][j].is_interior) continue;
          
          get_grid_neighbors_4(i, j, ni, nj);
          const std::size_t idxi = idx_map.at(std::make_pair(i, j));

          for (std::size_t k = 0; k < 4; ++k) { 
            const Size_pair pair = std::make_pair(ni[k], nj[k]);
            if (idx_map.find(pair) != idx_map.end()) {

              const std::size_t idxj = idx_map.at(pair);
              edges.push_back(std::make_pair(idxi, idxj));
              
              const double edge_weight = create_edge_weight(
                i, j, ni[k], nj[k], image);
              edge_weights.push_back(edge_weight);

            } else {
              // boundary pixel
            }
          }
        }
      }
    }

    void get_grid_neighbors_4(
      const std::size_t i, const std::size_t j,
      std::vector<std::size_t>& ni, 
      std::vector<std::size_t>& nj) {

      ni.clear(); nj.clear();
      ni.resize(4); nj.resize(4);

      CGAL_assertion(i > 0 && j > 0);

      ni[0] = i - 1; nj[0] = j;
      ni[1] = i;     nj[1] = j + 1;
      ni[2] = i + 1; nj[2] = j;
      ni[3] = i;     nj[3] = j - 1;
    }

    double create_edge_weight(
      const std::size_t i1, const std::size_t j1,
      const std::size_t i2, const std::size_t j2,
      const Image& image) {

      double edge_weight = 1.0;
      return CGAL::to_double(m_beta) * edge_weight;
    }

    void set_cost_matrix(
      const Image& image,
      const std::map<Size_pair, std::size_t>& idx_map,
      std::vector< std::vector<double> >& cost_matrix) {

      CGAL_assertion(idx_map.size() > 0);

      cost_matrix.clear();
      cost_matrix.resize(m_num_labels + 1); // -1
      for (std::size_t i = 0; i < m_num_labels + 1; ++i) // -1
        cost_matrix[i].resize(idx_map.size());

      std::vector<double> probabilities;
      for (std::size_t i = 1; i < image.rows - 1; ++i) {
        for (std::size_t j = 1; j < image.cols - 1; ++j) {
          // if (!image.grid[i][j].is_interior) continue;

          const std::size_t pixel_idx = idx_map.at(std::make_pair(i, j));

          create_probabilities(i, j, image, probabilities);
          for (std::size_t k = 0; k < m_num_labels + 1; ++k) // -1
            cost_matrix[k][pixel_idx] = 
              get_cost(image, i, j, probabilities[k]);
        }
      }
      // save_cost_matrix(image, idx_map, cost_matrix);
    }

    void create_probabilities(
      const std::size_t i, const std::size_t j,
      const Image& image,
      std::vector<double>& probabilities) {

      probabilities.clear();
      probabilities.resize(m_num_labels + 1, 0.0); // -1
      std::vector<std::size_t> nums(m_num_labels + 1, 0); // -1

      std::vector<std::size_t> ni, nj;
      get_grid_neighbors_8(i, j, ni, nj);

      for (std::size_t k = 0; k < 8; ++k) {
        const std::size_t ii = ni[k];
        const std::size_t jj = nj[k];

        const auto& cell = image.grid[ii][jj];
        // if (!cell.is_interior) continue;

        const std::size_t label = get_label(cell.zr, cell.zg, cell.zb);
        probabilities[label] += FT(1);
        nums[label] += 1;
      }

      double sum = 0.0;
      for (std::size_t k = 0; k < m_num_labels + 1; ++k) { // -1
        if (nums[k] == 0) continue;
        probabilities[k] /= static_cast<double>(nums[k]);
        sum += probabilities[k];
      }

      if (sum == 0.0)
        return;

      CGAL_assertion(sum > 0.0); double final_sum = 0.0;
      for (std::size_t k = 0; k < m_num_labels + 1; ++k) { // -1
        probabilities[k] /= sum;
        final_sum += probabilities[k];
      }
      CGAL_assertion(CGAL::abs(1.0 - final_sum) < 0.00001);
    }

    void get_grid_neighbors_8(
      const std::size_t i, const std::size_t j,
      std::vector<std::size_t>& ni, 
      std::vector<std::size_t>& nj) {

      ni.clear(); nj.clear();
      ni.resize(8); nj.resize(8);

      CGAL_assertion(i > 0 && j > 0);

      ni[0] = i - 1; nj[0] = j - 1;
      ni[1] = i - 1; nj[1] = j;
      ni[2] = i - 1; nj[2] = j + 1;
      ni[3] = i;     nj[3] = j + 1;
      ni[4] = i + 1; nj[4] = j + 1;
      ni[5] = i + 1; nj[5] = j;
      ni[6] = i + 1; nj[6] = j - 1;
      ni[7] = i;     nj[7] = j - 1;
    }

    double get_cost(
      const Image& image,
      const std::size_t i, const std::size_t j,
      const double prob) {
      
      const double weight = get_weight(i, j, image);
      return (1.0 - prob) * weight;
    }

    double get_weight(
      const std::size_t i, const std::size_t j,
      const Image& image) {

      return 1.0;
    }

    void compute_graphcut(
      const std::vector<Size_pair>& edges,
      const std::vector<double>& edge_weights,
      const std::vector< std::vector<double> >& cost_matrix,
      std::vector<std::size_t>& labels) {

      std::cout << "Initial labels (size " << 
      labels.size() << ")" << std::endl;

      Alpha_expansion graphcut;
      graphcut(edges, edge_weights, cost_matrix, labels);

      std::cout << "Final labels (size " << 
      labels.size() << ")" << std::endl;
    }

    bool is_outer_boundary_pixel(
      const long i1, const long j1,
      const long i2, const long j2) {

      const auto& cell1 = m_image.grid[i1][j1];
      const auto& cell2 = m_image.grid[i2][j2];

      return ( cell1.is_interior && !cell2.is_interior );
    }

    void save_cluster(const std::string name) {
      
      std::vector<Point_3> points;
      points.reserve(m_cluster.size());
      for (const auto& item: m_cluster)
        points.push_back(item.input_point);
      const Color color(0, 0, 0);
      m_saver.export_points(points, color, name);
    }

    void save_grid(const std::string name) {

      std::vector<Point_3> tmp;
      std::vector< std::vector<Point_3> > points;
      points.reserve(m_grid.size());

      for (const auto& pair : m_grid) {
        tmp.clear();
        for (const std::size_t idx : pair.second)
          tmp.push_back(m_cluster[idx].final_point);
        points.push_back(tmp);
      }
      m_saver.clear();
      m_saver.export_points(points, name);
    }

    void save_image(
      const std::string name,
      const Image& image) {
      
      OpenCVImage cvimage(
        image.rows * m_pixels_per_cell, 
        image.cols * m_pixels_per_cell, 
        CV_8UC3, cv::Scalar(255, 255, 255));

      for (std::size_t i = 0; i < image.rows; ++i) {
        for (std::size_t j = 0; j < image.cols; ++j) {
          
          const uchar zr = saturate_z(image.grid[i][j].zr);
          const uchar zg = saturate_z(image.grid[i][j].zg);
          const uchar zb = saturate_z(image.grid[i][j].zb);
          create_pixel(i, j, zr, zg, zb, cvimage);
        }
      }
      save_opencv_image(name, cvimage);
    }

    uchar saturate_z(const FT val) {
      const float z = static_cast<float>(val);
      return cv::saturate_cast<uchar>(z);
    }

    void create_pixel(
      const std::size_t i, const std::size_t j, 
      const uchar zr, const uchar zg, const uchar zb, 
      OpenCVImage& image) {

      const std::size_t il = i * m_pixels_per_cell;
      const std::size_t jl = j * m_pixels_per_cell;
      for (std::size_t ii = il; ii < il + m_pixels_per_cell; ++ii) {
        for (std::size_t jj = jl; jj < jl + m_pixels_per_cell; ++jj) {
          cv::Vec3b& bgr = image.at<cv::Vec3b>(ii, jj);
          bgr[0] = zb;
          bgr[1] = zg;
          bgr[2] = zr;
        }
      }
    }

    void save_opencv_image(
      const std::string& name,
      const OpenCVImage& image) {
      imwrite(name, image);
    }

    void save_point_cloud(
      const std::string name,
      const Image& image) {

      std::vector<Pixel> point_cloud;
      create_point_cloud(image, point_cloud);

      std::vector<Point_3> points;
      points.reserve(point_cloud.size());
      for (const auto& pixel : point_cloud) {
        if (!pixel.is_interior) continue;
        points.push_back(pixel.point);
      }
        
      const Color color(0, 0, 0);
      m_saver.export_points(points, color, name);
    }

    void save_cost_matrix(
      const Image& image,
      const std::map<Size_pair, std::size_t>& idx_map,
      const std::vector< std::vector<double> >& cost_matrix) {

      std::vector<Image> images(m_num_labels + 1); // -1
      for (std::size_t k = 0; k < m_num_labels + 1; ++k) // -1
        images[k].resize(image.rows, image.cols);

      for (std::size_t i = 1; i < image.rows - 1; ++i) {
        for (std::size_t j = 1; j < image.cols - 1; ++j) {
          // if (!image.grid[i][j].is_interior) continue;
          
          const std::size_t pixel_idx = idx_map.at(std::make_pair(i, j));
          for (std::size_t k = 0; k < m_num_labels + 1; ++k) { // -1
            const double prob = get_probability(cost_matrix[k][pixel_idx]);

            // Default color.
            images[k].grid[i][j].zr = FT(125);
            images[k].grid[i][j].zg = FT(0);
            images[k].grid[i][j].zb = FT(0);

            const FT z = FT(prob) * FT(255);
            images[k].create_pixel(i, j, 0, true, z, z, z);
          }
        }
      }

      for (std::size_t k = 0; k < m_num_labels + 1; ++k) { // -1
        const std::string name = "/Users/monet/Documents/lod/logs/buildings/tmp/" 
        + std::to_string(k) + "-image-probs.jpg";
        save_image(name, images[k]);
      }
    }

    double get_probability(const double cost) {
      return 1.0 - cost;
    }

    void save_regular_points(
      const std::vector< std::pair<Point_2, bool> >& input,
      const std::string name) {

      std::vector<Point_3> points;
      points.reserve(input.size());
      for (const auto& p : input)
        points.push_back(Point_3(p.first.x(), p.first.y(), FT(0)));
        
      const Color color(0, 0, 0);
      m_saver.export_points(points, color, name);
    }
  };

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_INTERNAL_GENERIC_SIMPLIFIER_H
