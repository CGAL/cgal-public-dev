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

#ifndef CGAL_LEVELS_OF_DETAIL_INTERNAL_CLOUD_TO_IMAGE_CONVERTER_H
#define CGAL_LEVELS_OF_DETAIL_INTERNAL_CLOUD_TO_IMAGE_CONVERTER_H

// STL includes.
#include <vector>
#include <utility>
#include <memory>
#include <stdio.h>

// CGAL includes.
#include <CGAL/assertions.h>
#include <CGAL/barycenter.h>
#include <CGAL/property_map.h>

#define CGAL_DO_NOT_USE_BOYKOV_KOLMOGOROV_MAXFLOW_SOFTWARE
#include <CGAL/internal/Surface_mesh_segmentation/Alpha_expansion_graph_cut.h>

// LOD includes.
#include <CGAL/Levels_of_detail/enum.h>

// Internal includes.
#include <CGAL/Levels_of_detail/internal/utils.h>
#include <CGAL/Levels_of_detail/internal/struct.h>

// Spatial search.
#include <CGAL/Levels_of_detail/internal/Spatial_search/K_neighbor_query.h>

// Shape detection.
#include <CGAL/Levels_of_detail/internal/Shape_detection/Estimate_normals_3.h>

// Buildings.
#include <CGAL/Levels_of_detail/internal/Buildings/Building_roofs_creator.h>

// Testing.
#include "../../../../../test/Levels_of_detail/include/Saver.h"

// Simplification.
#include <CGAL/Levels_of_detail/internal/Simplification/Alpha_shapes_filtering_2.h>

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
  typename PointMap2,
  typename PointMap3>
  class Cloud_to_image_converter {

  public:
    using Traits = GeomTraits;
    using Point_map_2 = PointMap2;
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
    using Double_pair = std::pair<double, std::size_t>;

    using K_neighbor_query = 
    internal::K_neighbor_query<Traits, Indices, Point_map_3>;
    
    using Normal_estimator_3 = 
    internal::Estimate_normals_3<Traits, Indices, Point_map_3, K_neighbor_query>;

    using Building_roofs_creator = 
    internal::Building_roofs_creator<Traits, Point_map_3>;

    using Identity_map_2 = 
    CGAL::Identity_property_map<Point_2>;

    using Alpha_shapes_filtering_2 = 
    internal::Alpha_shapes_filtering_2<Traits>;

    using Alpha_expansion = CGAL::internal::Alpha_expansion_graph_cut_boost;

    struct Cluster_item {
      Cluster_item(
        const Point_3 _point, 
        const std::size_t _idx) :
      input_point(_point),
      idx(_idx) { }
      
      Point_3 input_point;
      Point_3 final_point;
      std::size_t idx;
    };

    using Cell_id = std::pair<long, long>;
    using Cell_data = std::vector<std::size_t>;
    using Grid = std::map<Cell_id, Cell_data>;

    struct Image_cell {
      
      Image_cell() : 
      roof_idx(-1),
      zr(FT(125)), zg(FT(0)), zb(FT(0)),
      is_interior(false) { }

      std::size_t roof_idx;
      FT zr, zg, zb;
      bool is_interior;
    };

    struct Image {

      Image() : rows(0), cols(0) { }
      Image(const long _rows, const long _cols) : 
      rows(_rows), cols(_cols) { 
        resize(_rows, _cols);
      }

      void resize(const long _rows, const long _cols) {
        rows = _rows; cols = _cols;
        grid.resize(_rows);
        for (auto& pixels : grid)
          pixels.resize(_cols);
      }

      void create_pixel(
        const long i, const long j,
        const std::size_t roof_idx,
        const FT zr, const FT zg, const FT zb) {

        auto& pixel = grid[i][j];

        pixel.roof_idx = roof_idx;
        pixel.zr = zr;
        pixel.zg = zg;
        pixel.zb = zb;
        pixel.is_interior = true;
      }

      void recolor_pixel(
        const long i, const long j,
        const FT zr, const FT zg, const FT zb) {
        
        auto& pixel = grid[i][j];
        pixel.zr = zr;
        pixel.zg = zg;
        pixel.zb = zb; 
        pixel.is_interior = true;
      }

      long rows, cols;
      std::vector< std::vector<Image_cell> > grid;
    };

    using OpenCVImage = cv::Mat;

    Cloud_to_image_converter(
      const Indices& input_range,
      const Point_map_2 point_map_2,
      const Point_map_3 point_map_3,
      const FT grid_cell_width_2 = 0.25,
      const FT region_growing_scale_3 = 0.8,
      const FT region_growing_noise_level_3 = 0.5,
      const FT region_growing_angle_3 = 25.0,
      const FT region_growing_min_area_3 = 4.0,
      const FT region_growing_distance_to_line_3 = 0.25,
      const FT alpha_shape_size_2 = 0.5,
      const FT beta = 1.0) :
    m_input_range(input_range),
    m_point_map_2(point_map_2),
    m_point_map_3(point_map_3),
    m_val_min(+internal::max_value<FT>()),
    m_val_max(-internal::max_value<FT>()),
    m_rows_min(+internal::max_value<long>()),
    m_rows_max(-internal::max_value<long>()),
    m_cols_min(+internal::max_value<long>()),
    m_cols_max(-internal::max_value<long>()),
    m_samples_per_face(20), // should be in [0,100]
    m_pixels_per_cell(9), // should be an odd number, change later
    m_lsd_scale(FT(9) / FT(10)), // should be in [0,1)
    // grid parameters:
    m_grid_cell_width_2(grid_cell_width_2),
    // region growing parameters:
    m_region_growing_scale_3(region_growing_scale_3),
    m_region_growing_noise_level_3(region_growing_noise_level_3),
    m_region_growing_angle_3(region_growing_angle_3),
    m_region_growing_min_area_3(region_growing_min_area_3),
    m_region_growing_distance_to_line_3(region_growing_distance_to_line_3),
    m_alpha_shape_size_2(alpha_shape_size_2),
    m_beta(beta),
    m_num_labels(0)
    { }

    void convert() {

      create_cluster();
      transform_cluster();
      create_grid();
      create_image();
      // apply_lsd();
    }

  private:
    const Indices& m_input_range;
    const Point_map_2 m_point_map_2;
    const Point_map_3 m_point_map_3;

    Indices m_clean_input;
    std::vector<Indices> m_roof_regions;
    std::vector<Points_3> m_roofs;
    std::map<std::size_t, FT> m_height_map;
    std::vector<Cluster_item> m_cluster;
    FT m_val_min, m_val_max;
    
    Grid m_grid;
    long m_rows_min, m_rows_max;
    long m_cols_min, m_cols_max;
    
    const std::size_t m_samples_per_face;
    const long m_pixels_per_cell;
    const FT m_lsd_scale;

    const FT m_grid_cell_width_2;
    const FT m_region_growing_scale_3;
    const FT m_region_growing_noise_level_3;
    const FT m_region_growing_angle_3;
    const FT m_region_growing_min_area_3;
    const FT m_region_growing_distance_to_line_3;
    const FT m_alpha_shape_size_2;
    const FT m_beta;

    std::size_t m_num_labels;

    Saver m_saver;

    void create_cluster() {
      find_roof_regions();
      create_roofs();
      add_roof_points_to_cluster();
      save_cluster("/Users/monet/Documents/lod/logs/buildings/cluster");
    }

    void find_roof_regions() {

      const Building_roofs_creator creator(m_point_map_3);

      CGAL_assertion(m_input_range.size() >= 3);
      creator.create_cluster(
        m_input_range,
        m_region_growing_scale_3,
        m_region_growing_angle_3,
        m_clean_input);

      CGAL_assertion(m_clean_input.size() >= 3);
      creator.create_roof_regions(
        m_clean_input,
        m_region_growing_scale_3,
        m_region_growing_noise_level_3,
        m_region_growing_angle_3,
        m_region_growing_min_area_3,
        m_region_growing_distance_to_line_3,
        m_alpha_shape_size_2,
        m_roof_regions);
    }

    void create_roofs() {

      m_roofs.clear();
      m_roofs.reserve(m_roof_regions.size());
      m_height_map.clear();

      Points_3 roof; Plane_3 plane;
      for (std::size_t i = 0; i < m_roof_regions.size(); ++i) {
        const auto& roof_region = m_roof_regions[i];
        roof.clear();
  
        internal::plane_from_points_3(
        m_clean_input, m_point_map_3, roof_region, plane);
        internal::project_on_plane_3(
        m_clean_input, m_point_map_3, roof_region, plane, roof);
        sample_roof_region(plane, roof);
        m_roofs.push_back(roof);
        
        const FT height = get_roof_height(roof);
        m_height_map[i] = height;
      }
    }

    FT get_roof_height(const Points_3& roof) const {

      FT height = -internal::max_value<FT>();
      for (const auto& p : roof)
        height = CGAL::max(height, p.z());
      return height;
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

    void add_roof_points_to_cluster() {

      CGAL_assertion(m_roofs.size() >= 1);
      m_cluster.clear();
      for (std::size_t i = 0; i < m_roofs.size(); ++i) {
        const auto& roof = m_roofs[i];
        for (const auto& point : roof) {
          m_cluster.push_back(Cluster_item(point, i));
          m_val_min = CGAL::min(point.z(), m_val_min);
          m_val_max = CGAL::max(point.z(), m_val_max);
        }
      }
    }

    void transform_cluster() {

      std::vector<Point_2> points;
      points.reserve(m_cluster.size());

      for (const auto& item : m_cluster)
        points.push_back(internal::point_2_from_point_3(item.input_point));

      Vector_2 dir;
      internal::estimate_direction_2(points, dir);
      const Vector_2 y_dir = Vector_2(FT(0), FT(1));

      FT angle_2d;
			internal::compute_angle_2(dir, y_dir, angle_2d);

      Point_2 b;
      internal::compute_barycenter_2(points, b);
      
      for (Point_2& p : points)
        internal::rotate_point_2(angle_2d, b, p);

      CGAL::Identity_property_map<Point_2> pmap;
      std::vector<Point_2> bbox;
      internal::bounding_box_2(points, pmap, bbox);

      const Point_2& tr = bbox[0];
      for (Point_2& p : points) {
        const FT x = p.x() - tr.x();
        const FT y = p.y() - tr.y();
        p = Point_2(x, y);
      }

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
      save_grid("/Users/monet/Documents/lod/logs/buildings/grid");
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

    void create_image() {

      const long rowsdiff = m_rows_max - m_rows_min;
      const long colsdiff = m_cols_max - m_cols_min;
      const long rows = (rowsdiff + 3); // +2 = +1 (diff pixel) +2 (boundary pixels)
      const long cols = (colsdiff + 3); // +2 = +1 (diff pixel) +2 (boundary pixels)

      std::cout << "Resolution (original): " << cols << "x" << rows << std::endl;
      std::cout << "Cols: " << colsdiff << " Rows: " << rowsdiff << std::endl;
      std::cout << "Val min: " << m_val_min << " Val max: " << m_val_max << std::endl;

      Image image(rows, cols);

      initialize_image(image);
      save_image("/Users/monet/Documents/lod/logs/buildings/image-origin.jpg", image);
      
      inpaint_image_opencv(image);
      save_image("/Users/monet/Documents/lod/logs/buildings/image-paints.jpg", image);

      apply_graphcut(image);
      save_image("/Users/monet/Documents/lod/logs/buildings/image-gcuted.jpg", image);

      // Interior points.
      // save_image("/Users/monet/Documents/lod/logs/buildings/image-green.jpg", image, true);
    }

    void initialize_image(
      Image& image) const {

      init_image(image);
    }

    void init_image(
      Image& image) const {

      long numcells = 0;
      for (long i = 1; i < image.rows - 1; ++i) {
        for (long j = 1; j < image.cols - 1; ++j) {

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

    long get_id_x(const long j) const {
      return m_cols_min + j;
    }

    long get_id_y(const long i) const {
      return m_rows_max - i;
    }

    void initialize_pixel(
      const long i, const long j,
      const Cell_data& indices,
      Image& image) const {

      initialize_pixel_max_height(
        i, j, indices, image);
    }

    void initialize_pixel_max_height(
      const long i, const long j,
      const Cell_data& indices,
      Image& image) const {
    
      const auto& pair = max_z_height(indices);
      init_pixel(i, j, pair.first, pair.second, image);
    }

    std::pair<std::size_t, FT> max_z_height(const Cell_data& indices) const {

      std::map<std::size_t, std::size_t> vals;
      for (const std::size_t idx : indices)
        vals[m_cluster[idx].idx] = FT(0);
      for (const std::size_t idx : indices)
        vals[m_cluster[idx].idx] += FT(1);

      FT maxv = -FT(1);
      std::size_t final_idx = 0;
      for (const auto& pair : vals) {
        if (pair.second > maxv) {
          final_idx = pair.first;
          maxv = pair.second;
        }
      }

      const FT height = m_height_map.at(final_idx);
      return std::make_pair(final_idx, height);
    }

    void initialize_pixel_naive(
      const long i, const long j,
      const Cell_data& indices,
      Image& image) const {
      
      const FT val = average_z_height(indices);
      init_pixel(i, j, 0, val, image);
    }

    void init_pixel(
      const long i, const long j,
      const std::size_t roof_idx,
      const FT val,
      Image& image) const {
    
      const FT nor = normalize_z(val);
      image.create_pixel(i, j, roof_idx, nor, nor, nor);
    }

    FT average_z_height(const Cell_data& indices) const {
      CGAL_assertion(indices.size() > 0);
      FT val = FT(0);
      for (const std::size_t idx : indices)
        val += m_cluster[idx].final_point.z();
      val /= indices.size();
      return val;
    }

    FT normalize_z(const FT val) const {
      const FT z = (val - m_val_min) / (m_val_max - m_val_min);
      return z * FT(255);
    }

    void inpaint_image_opencv(
      Image& image) {
      
      OpenCVImage input(
        image.rows, 
        image.cols, 
        CV_8UC3, cv::Scalar(0, 0, 125));

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

      Mat inpainted;
      inpaint(input, mask, inpainted, 0, cv::INPAINT_TELEA);

      Image colored(image.rows, image.cols);
      for (std::size_t i = 1; i < colored.rows - 1; ++i) {
        for (std::size_t j = 1; j < colored.cols - 1; ++j) {
          const cv::Vec3b& bgr = inpainted.at<cv::Vec3b>(i, j);
          const std::size_t roof_idx = image.grid[i][j].roof_idx;

          const FT zr = FT(static_cast<std::size_t>(bgr[2]));
          const FT zg = FT(static_cast<std::size_t>(bgr[1]));
          const FT zb = FT(static_cast<std::size_t>(bgr[0]));

          colored.create_pixel(i, j, roof_idx, zr, zg, zb);
        }
      }
      image = colored;
    }

    void apply_graphcut(
      Image& image) {

      std::map<std::size_t, std::size_t> label_map, inv_label_map;
      set_label_map(image, label_map, inv_label_map);

      std::map<Size_pair, std::size_t> idx_map;
      set_idx_map(image, idx_map);

      std::vector<std::size_t> labels;
      set_initial_labels(image, idx_map, label_map, inv_label_map, labels);
      apply_new_labels(idx_map, inv_label_map, labels, image);
      save_image("/Users/monet/Documents/lod/logs/buildings/image-labels.jpg", image);

      std::vector<Size_pair> edges;
      std::vector<double> edge_weights;
      set_graphcut_edges(image, idx_map, edges, edge_weights);
      
      std::vector< std::vector<double> > cost_matrix;
      set_cost_matrix(image, idx_map, label_map, cost_matrix);

      compute_graphcut(edges, edge_weights, cost_matrix, labels);
      apply_new_labels(idx_map, inv_label_map, labels, image);
    }

    void set_label_map(
      const Image& image,
      std::map<std::size_t, std::size_t>& label_map,
      std::map<std::size_t, std::size_t>& inv_label_map) {
      
      label_map.clear();
      for (std::size_t i = 1; i < image.rows - 1; ++i) {
        for (std::size_t j = 1; j < image.cols - 1; ++j) {
          const std::size_t key = get_key(image.grid[i][j].zr);
          if (image.grid[i][j].roof_idx != std::size_t(-1))
            label_map[key] = key;
        }
      }
      m_num_labels = label_map.size();
      std::cout << "Num labels: " << m_num_labels << std::endl;
      
      std::size_t count = 0;
      for (auto& pair : label_map) {
        inv_label_map[count] = pair.first;
        pair.second = count;
        ++count;
      }
    }

    std::size_t get_key(const FT val) const {
      return static_cast<std::size_t>(CGAL::to_double(val)) / 10;
    }

    void set_idx_map(
      const Image& image,
      std::map<Size_pair, std::size_t>& idx_map) const {

      idx_map.clear();
      std::size_t pixel_idx = 0;
      for (std::size_t i = 1; i < image.rows - 1; ++i) {
        for (std::size_t j = 1; j < image.cols - 1; ++j) {
          idx_map[std::make_pair(i, j)] = pixel_idx;
          ++pixel_idx;
        }
      }
    }

    void set_graphcut_edges(
      const Image& image,
      const std::map<Size_pair, std::size_t>& idx_map,
      std::vector<Size_pair>& edges,
      std::vector<double>& edge_weights) const {

      edges.clear();
      edge_weights.clear();
      std::vector<std::size_t> ni, nj;
      for (std::size_t i = 1; i < image.rows - 1; ++i) {
        for (std::size_t j = 1; j < image.cols - 1; ++j) {
          
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
      std::vector<std::size_t>& nj) const {

      ni.clear(); nj.clear();
      ni.resize(4); nj.resize(4);

      CGAL_assertion(i > 0 && j > 0);

      ni[0] = i - 1; nj[0] = j;
      ni[1] = i;     nj[1] = j + 1;
      ni[2] = i + 1; nj[2] = j;
      ni[3] = i;     nj[3] = j - 1;
    }

    void get_grid_neighbors_8(
      const std::size_t i, const std::size_t j,
      std::vector<std::size_t>& ni, 
      std::vector<std::size_t>& nj) const {

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

    double create_edge_weight(
      const std::size_t i1, const std::size_t j1,
      const std::size_t i2, const std::size_t j2,
      const Image& image) const {

      // const FT val1 = image.grid[i1][j1].zr;
      // const FT val2 = image.grid[i2][j2].zr;
      // if (val1 != val2) edge_weight *= 0.0;

      double edge_weight = 1.0;
      return CGAL::to_double(m_beta) * edge_weight;
    }

    void set_cost_matrix(
      const Image& image,
      const std::map<Size_pair, std::size_t>& idx_map,
      const std::map<std::size_t, std::size_t>& label_map,
      std::vector< std::vector<double> >& cost_matrix) {

      CGAL_assertion(idx_map.size() > 0);
      CGAL_assertion(m_num_labels >= 2);

      cost_matrix.clear();
      cost_matrix.resize(m_num_labels);
      for (std::size_t i = 0; i < m_num_labels; ++i)
        cost_matrix[i].resize(idx_map.size());

      std::vector<double> probabilities;
      for (std::size_t i = 1; i < image.rows - 1; ++i) {
        for (std::size_t j = 1; j < image.cols - 1; ++j) {
          const std::size_t pixel_idx = idx_map.at(std::make_pair(i, j));

          create_probabilities(i, j, image, label_map, probabilities);
          for (std::size_t k = 0; k < m_num_labels; ++k)
            cost_matrix[k][pixel_idx] = get_cost(probabilities[k]);
        }
      }
      // save_cost_matrix(image, idx_map, cost_matrix);
    }

    double get_cost(const double prob) const {
      return 1.0 - prob;
    }

    double get_probability(const double cost) const {
      return 1.0 - cost;
    }

    void save_cost_matrix(
      const Image& image,
      const std::map<Size_pair, std::size_t>& idx_map,
      const std::vector< std::vector<double> >& cost_matrix) {

      std::vector<Image> images(m_num_labels);
      for (std::size_t k = 0; k < m_num_labels; ++k)
        images[k].resize(image.rows, image.cols);

      for (std::size_t i = 1; i < image.rows - 1; ++i) {
        for (std::size_t j = 1; j < image.cols - 1; ++j) {
          
          const std::size_t pixel_idx = idx_map.at(std::make_pair(i, j));
          for (std::size_t k = 0; k < m_num_labels; ++k) {
            const double prob = get_probability(cost_matrix[k][pixel_idx]);
            
            const FT z = FT(prob) * FT(255);
            images[k].create_pixel(i, j, 0, z, z, z);
          }
        }
      }

      for (std::size_t k = 0; k < m_num_labels; ++k) {
        const std::string name = "/Users/monet/Documents/lod/logs/buildings/images/" 
        + std::to_string(k) + "-image-probs.jpg";
        save_image(name, images[k]);
      }
    }

    void create_probabilities(
      const std::size_t i, const std::size_t j,
      const Image& image,
      const std::map<std::size_t, std::size_t>& label_map,
      std::vector<double>& probabilities) const {

      probabilities.clear();
      probabilities.resize(m_num_labels, 0.0);
      std::vector<std::size_t> nums(m_num_labels, 0);

      std::vector<std::size_t> ni, nj;
      get_grid_neighbors_8(i, j, ni, nj);

      for (std::size_t k = 0; k < 8; ++k) {
        const std::size_t ii = ni[k];
        const std::size_t jj = nj[k];

        const auto& cell = image.grid[ii][jj];
        const FT value = cell.zr;
        const std::size_t key = get_key(value);
        const std::size_t idx = label_map.at(key);
        const double prob = CGAL::to_double(value / FT(255));
        CGAL_assertion(idx >= 0 && idx < m_num_labels);
        probabilities[idx] += prob;
        nums[idx] += 1;
      }

      double sum = 0.0;
      for (std::size_t k = 0; k < m_num_labels; ++k) {
        if (nums[k] == 0) continue;
        probabilities[k] /= static_cast<double>(nums[k]);
        sum += probabilities[k];
      }

      CGAL_assertion(sum > 0.0); double final_sum = 0.0;
      for (std::size_t k = 0; k < m_num_labels; ++k) {
        probabilities[k] /= sum;
        final_sum += probabilities[k];
      }
      CGAL_assertion(CGAL::abs(1.0 - final_sum) < 0.00001);
    }

    void set_initial_labels(
      const Image& image,
      const std::map<Size_pair, std::size_t>& idx_map,
      const std::map<std::size_t, std::size_t>& label_map,
      const std::map<std::size_t, std::size_t>& inv_label_map,
      std::vector<std::size_t>& labels) {

      labels.clear();
      labels.resize(idx_map.size());
      for (std::size_t i = 1; i < image.rows - 1; ++i) {
        for (std::size_t j = 1; j < image.cols - 1; ++j) {
          
          std::size_t key = get_key(image.grid[i][j].zr);
          update_key(label_map, key);
          
          const std::size_t pixel_idx = idx_map.at(std::make_pair(i, j));
          labels[pixel_idx] = label_map.at(key);
        }
      }
    }

    void update_key(
      const std::map<std::size_t, std::size_t>& label_map,
      std::size_t& label) const {

      std::size_t mindiff = internal::max_value<std::size_t>();
      std::size_t closest = label;

      for (const auto& pair : label_map) {
        const std::size_t key = pair.first;
        
        std::size_t diff;
        if (key >= label) diff = key - label;
        else diff = label - key;

        if (diff < mindiff) {
          closest = key;
          mindiff = diff;
        }
      }
      label = closest;
    }

    void compute_graphcut(
      const std::vector<Size_pair>& edges,
      const std::vector<double>& edge_weights,
      const std::vector< std::vector<double> >& probability_matrix,
      std::vector<std::size_t>& labels) const {

      std::cout << "Initial labels (size " << 
      labels.size() << ")" << std::endl;

      Alpha_expansion graphcut;
      graphcut(edges, edge_weights, probability_matrix, labels);

      std::cout << "Final labels (size " << 
      labels.size() << ")" << std::endl;
    }

    void apply_new_labels(
      const std::map<Size_pair, std::size_t>& idx_map,
      const std::map<std::size_t, std::size_t>& inv_label_map,
      const std::vector<std::size_t>& labels,
      Image& image) {

      Image labeled(image.rows, image.cols);
      for (std::size_t i = 1; i < image.rows - 1; ++i) {
        for (std::size_t j = 1; j < image.cols - 1; ++j) {
          
          const std::size_t pixel_idx = idx_map.at(std::make_pair(i, j));
          const std::size_t z = inv_label_map.at(labels[pixel_idx]);

          const FT zr = FT(z * 10);
          const FT zg = FT(z * 10);
          const FT zb = FT(z * 10);
          
          const std::size_t roof_idx = image.grid[i][j].roof_idx;
          labeled.create_pixel(i, j, roof_idx, zr, zg, zb);
        }
      }
      image = labeled;
    }

    void apply_lsd() {
      
      OpenCVImage image;
      read_gray_scale_opencv_image(
        "/Users/monet/Documents/lod/logs/buildings/image-gcuted.jpg", image);
      save_opencv_image(
        "/Users/monet/Documents/lod/logs/buildings/image-gscale.jpg", image);

      vector<cv::Vec4f> lines_std;
      cv::Ptr<cv::LineSegmentDetector> ls = createLineSegmentDetector(
        cv::LSD_REFINE_STD, m_lsd_scale);
      ls->detect(image, lines_std);

      save_opencv_lines(
        "/Users/monet/Documents/lod/logs/buildings/result-lsd.jpg", 
        image, ls, lines_std);
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

    void read_gray_scale_opencv_image(
      const std::string name,
      OpenCVImage& image) {

      image = imread(name, cv::IMREAD_GRAYSCALE);
    }

    void save_opencv_image(
      const std::string& name,
      const OpenCVImage& image) {

      imwrite(name, image);
    }

    void save_image(
      const std::string name,
      const Image& image,
      const bool interior = false) {
      
      OpenCVImage cvimage(
        image.rows * m_pixels_per_cell, 
        image.cols * m_pixels_per_cell, 
        CV_8UC3, cv::Scalar(0, 0, 125));
      // std::cout << "Resolution (opencv): " << 
      // cvimage.cols << "x" << cvimage.rows << std::endl;

      for (long i = 0; i < image.rows; ++i) {
        for (long j = 0; j < image.cols; ++j) {
          
          if (!interior) {
            const uchar zr = saturate_z(image.grid[i][j].zr);
            const uchar zg = saturate_z(image.grid[i][j].zg);
            const uchar zb = saturate_z(image.grid[i][j].zb);
            create_pixel(i, j, zr, zg, zb, cvimage);

          } else {

            if (image.grid[i][j].is_interior) {
              const uchar zr = saturate_z(0.0);
              const uchar zg = saturate_z(125.0);
              const uchar zb = saturate_z(0.0);
              create_pixel(i, j, zr, zg, zb, cvimage);

            } else {
              const uchar zr = saturate_z(125.0);
              const uchar zg = saturate_z(0.0);
              const uchar zb = saturate_z(0.0);
              create_pixel(i, j, zr, zg, zb, cvimage);
            }
          }
        }
      }
      save_opencv_image(name, cvimage);
    }

    void create_pixel(
      const long i, const long j, 
      const uchar zr, const uchar zg, const uchar zb, 
      OpenCVImage& image) const {

      const long il = i * m_pixels_per_cell;
      const long jl = j * m_pixels_per_cell;
      for (long ii = il; ii < il + m_pixels_per_cell; ++ii) {
        for (long jj = jl; jj < jl + m_pixels_per_cell; ++jj) {
          cv::Vec3b& bgr = image.at<cv::Vec3b>(ii, jj);
          bgr[0] = zb;
          bgr[1] = zg;
          bgr[2] = zr;
        }
      }
    }

    uchar saturate_z(const FT val) const {
      const float z = static_cast<float>(val);
      return cv::saturate_cast<uchar>(z);
    }

    void save_opencv_lines(
      const std::string name,
      const OpenCVImage& image,
      const cv::Ptr<cv::LineSegmentDetector>& ls,
      const vector<cv::Vec4f>& lines_std) {

      OpenCVImage drawn_lines(image); 
      ls->drawSegments(drawn_lines, lines_std);
      save_opencv_image(name, drawn_lines);
    }
  };

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_INTERNAL_CLOUD_TO_IMAGE_CONVERTER_H
