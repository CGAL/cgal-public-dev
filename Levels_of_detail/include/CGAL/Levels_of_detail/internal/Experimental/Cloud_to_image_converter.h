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
      tagged(false) { }

      std::size_t roof_idx;
      FT zr, zg, zb;
      bool tagged;
    };

    struct Image {

      Image() : rows(0), cols(0) { }
      Image(const long _rows, const long _cols) : 
      rows(_rows), cols(_cols) { 
        grid.resize(rows);
        for (auto& pixels : grid)
          pixels.resize(cols);
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
        pixel.tagged = true;
      }

      void recolor_pixel(
        const long i, const long j,
        const FT zr, const FT zg, const FT zb) {
        
        auto& pixel = grid[i][j];
        pixel.zr = zr;
        pixel.zg = zg;
        pixel.zb = zb; 
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
      const FT alpha_shape_size_2 = 0.5) :
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
    m_alpha_shape_size_2(alpha_shape_size_2)
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
      const long rows = (rowsdiff + 1);
      const long cols = (colsdiff + 1);

      std::cout << "Resolution (original): " << cols << "x" << rows << std::endl;
      std::cout << "Cols: " << colsdiff << " Rows: " << rowsdiff << std::endl;
      std::cout << "Val min: " << m_val_min << " Val max: " << m_val_max << std::endl;

      Image image(rows, cols);
      initialize_image(image);
      save_image("/Users/monet/Documents/lod/logs/buildings/image-origin.jpg", image);
      interpolate_values(image);
      save_image("/Users/monet/Documents/lod/logs/buildings/image-interp.jpg", image);
    }

    void initialize_image(
      Image& image) const {

      long numcells = 0;
      for (long i = 0; i < image.rows; ++i) {
        for (long j = 0; j < image.cols; ++j) {

          const long id_x = get_id_x(j);
          const long id_y = get_id_y(i);
          
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

      // initialize_pixel_naive(
      //   i, j, indices, image);
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
      
      const FT val  = average_z_height(indices);
      init_pixel(i, j, 0, val, image);
    }

    void init_pixel(
      const long i, const long j,
      const std::size_t roof_idx,
      const FT val,
      Image& image) const {
    
      const FT nor  = normalize_z(val);
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

    void interpolate_values(
      Image& image) const {

      close_one_pixel_holes(image);
      guess_values(image);
    }

    void close_one_pixel_holes(
      Image& image) const {

      std::vector<long> ni, nj;
      for (long i = 0; i < image.rows; ++i) {
        for (long j = 0; j < image.cols; ++j) {
          get_grid_neighbors(i, j, ni, nj);

          const FT zr_or = image.grid[i][j].zr;
          const FT zg_or = image.grid[i][j].zg;
          const FT zb_or = image.grid[i][j].zb;

          FT zr_ref = FT(-1), zg_ref = FT(-1), zb_ref = FT(-1);
          find_reference_color(image, ni, nj, zr_ref, zg_ref, zb_ref);
          CGAL_assertion(zr_ref != FT(-1) && zg_ref != FT(-1) && zb_ref != FT(-1));

          bool found = true;
          for (std::size_t k = 0; k < 8; ++k) {
            if (is_exterior_pixel(image.rows, image.cols, ni[k], nj[k]))
              continue;

            const FT zr = image.grid[ni[k]][nj[k]].zr;
            const FT zg = image.grid[ni[k]][nj[k]].zg;
            const FT zb = image.grid[ni[k]][nj[k]].zb;

            if (zr != zr_ref || zg != zg_ref || zb != zb_ref) {
              found = false;
              break;
            }
          }

          if (found) {
            if (zr_or != zr_ref || zg_or != zg_ref || zb_or != zb_ref)
              image.recolor_pixel(i, j, zr_ref, zg_ref, zb_ref);
          }
        }
      }
    }

    void find_reference_color(
      const Image& image,
      const std::vector<long>& ni,
      const std::vector<long>& nj,
      FT& zr, FT& zg, FT& zb) const {

      for (std::size_t k = 0; k < 8; ++k) {
        if (is_exterior_pixel(image.rows, image.cols, ni[k], nj[k]))
          continue;

        zr = image.grid[ni[k]][nj[k]].zr;
        zg = image.grid[ni[k]][nj[k]].zg;
        zb = image.grid[ni[k]][nj[k]].zb;
        return;
      }
    }

    void get_grid_neighbors(
      const long i, const long j,
      std::vector<long>& ni, 
      std::vector<long>& nj) const {

      ni.clear(); nj.clear();
      ni.resize(8); nj.resize(8);

      ni[0] = i - 1; nj[0] = j - 1;
      ni[1] = i - 1; nj[1] = j;
      ni[2] = i - 1; nj[2] = j + 1;
      ni[3] = i;     nj[3] = j + 1;
      ni[4] = i + 1; nj[4] = j + 1;
      ni[5] = i + 1; nj[5] = j;
      ni[6] = i + 1; nj[6] = j - 1;
      ni[7] = i;     nj[7] = j - 1;
    }

    bool is_exterior_pixel(
      const long rows, const long cols,
      const long i, const long j) const {

      return (i < 0 || i >= rows || j < 0 || j >= cols);
    }

    void guess_values(Image& image) const {

    }

    void apply_lsd() {
      
      OpenCVImage image;
      read_gray_scale_opencv_image(
        "/Users/monet/Documents/lod/logs/buildings/image-interp.jpg", image);
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
      const Image& image) {
      
      OpenCVImage cvimage(
        image.rows * m_pixels_per_cell, 
        image.cols * m_pixels_per_cell, 
        CV_8UC3, cv::Scalar(0, 0, 125));
      std::cout << "Resolution (opencv): " << 
      cvimage.cols << "x" << cvimage.rows << std::endl;

      for (long i = 0; i < image.rows; ++i) {
        for (long j = 0; j < image.cols; ++j) {
          const uchar zr = saturate_z(image.grid[i][j].zr);
          const uchar zg = saturate_z(image.grid[i][j].zg);
          const uchar zb = saturate_z(image.grid[i][j].zb);
          create_pixel(i, j, zr, zg, zb, cvimage);
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
