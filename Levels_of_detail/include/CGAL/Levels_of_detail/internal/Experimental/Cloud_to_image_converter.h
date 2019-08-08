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
#include <queue>
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
      rotate_cluster();
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
    std::vector<Cluster_item> m_cluster;
    FT m_val_min, m_val_max;
    
    Grid m_grid;
    long m_rows_min, m_rows_max;
    long m_cols_min, m_cols_max;
    
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

      Points_3 roof; Plane_3 plane;
      for (const auto& roof_region : m_roof_regions) {
        roof.clear();
  
        internal::plane_from_points_3(
        m_clean_input, m_point_map_3, roof_region, plane);
        internal::project_on_plane_3(
        m_clean_input, m_point_map_3, roof_region, plane, roof);
        sample_roof_region(plane, roof);
        m_roofs.push_back(roof);
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
      points.clear(); filtering.get_samples(sampling_2, 20, points);
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

    void rotate_cluster() {

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

      // Save results.
      std::vector<Point_3> tmp;
      std::vector< std::vector<Point_3> > points;
      points.reserve(m_grid.size());

      for (const auto& pair : m_grid) {
        tmp.clear();
        for (const std::size_t idx : pair.second)
          tmp.push_back(m_cluster[idx].final_point);
        points.push_back(tmp);
      }

      const std::string name =
      "/Users/monet/Documents/lod/logs/buildings/grid";
      m_saver.export_points(points, name);
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

      const long diff1 = m_rows_max - m_rows_min;
      const long diff2 = m_cols_max - m_cols_min;

      const long size = CGAL::max(diff1, diff2) + 1;
      const long rows = size * m_pixels_per_cell;
      const long cols = size * m_pixels_per_cell;

      cv::Mat image(rows, cols, CV_8UC3, cv::Scalar(0, 0, 125));
      long numcells = 0;
      std::queue<long> iis, jjs;

      // Create an initial image with red bottom.
      for (long i = 0; i < image.rows / m_pixels_per_cell; ++i) {
        for (long j = 0; j < image.cols / m_pixels_per_cell; ++j) {

          const long id_x = m_cols_min + i;
          const long id_y = m_rows_max - j;
          
          const Cell_id cell_id = std::make_pair(id_x, id_y);
          if (m_grid.find(cell_id) != m_grid.end()) {
            ++numcells;
            const auto& data = m_grid.at(cell_id);

            FT val = FT(0);
            for (const std::size_t idx : data)
              val += m_cluster[idx].final_point.z();
            val /= data.size();

            const FT norm = (val - m_val_min) / (m_val_max - m_val_min);
            const FT z = norm * FT(255);
            const float zf = static_cast<float>(z);
            const uchar finalz = cv::saturate_cast<uchar>(zf);

            const long il = i * m_pixels_per_cell;
            const long jl = j * m_pixels_per_cell;

            for (long ii = il; ii < il + m_pixels_per_cell; ++ii) {
              for (long jj = jl; jj < jl + m_pixels_per_cell; ++jj) {
                cv::Vec3b& bgr = image.at<cv::Vec3b>(ii, jj);

                bgr[0] = finalz;
                bgr[1] = finalz;
                bgr[2] = finalz;
              }
            }
          } else {
            iis.push(i);
            jjs.push(j);
          }
        }
      }

      // Interpolate through all unfilled entries, which are red.
      while (!iis.empty() && !jjs.empty()) {

        const long i = iis.front();
        const long j = jjs.front();
        iis.pop(); jjs.pop();

        // Get neighbors.
        std::vector<long> ni(8), nj(8);
        ni[0] = i - 1; nj[0] = j - 1;
        ni[1] = i; nj[1] = j - 1;
        ni[2] = i + 1; nj[2] = j - 1;
        ni[3] = i + 1; nj[3] = j;
        ni[4] = i + 1; nj[4] = j + 1;
        ni[5] = i; nj[5] = j + 1;
        ni[6] = i - 1; nj[6] = j + 1;
        ni[7] = i - 1; nj[7] = j;
        
        FT val = FT(0); FT numvals = FT(0);
        for (std::size_t k = 0; k < 8; ++k) {

          // Boundary index.
          if (ni[k] < 0 || ni[k] >= size || nj[k] < 0 || nj[k] >= size)
            continue;

          // All other indices.
          const long id_x = m_cols_min + ni[k];
          const long id_y = m_rows_max - nj[k];

          const Cell_id cell_id = std::make_pair(id_x, id_y);
          if (m_grid.find(cell_id) == m_grid.end())
            continue;
          
          // Add value.
          const auto& data = m_grid.at(cell_id);

          FT tval = FT(0);
          for (const std::size_t idx : data)
            tval += m_cluster[idx].final_point.z();
          tval /= data.size();
          val += tval;
          numvals += FT(1);
        }

        if (numvals < FT(1))
          continue;
        val /= numvals;

        const FT norm = (val - m_val_min) / (m_val_max - m_val_min);
        const FT z = norm * FT(255);
        const float zf = static_cast<float>(z);
        const uchar finalz = cv::saturate_cast<uchar>(zf);

        const long il = i * m_pixels_per_cell;
        const long jl = j * m_pixels_per_cell;

        for (long ii = il; ii < il + m_pixels_per_cell; ++ii) {
          for (long jj = jl; jj < jl + m_pixels_per_cell; ++jj) {
            cv::Vec3b& bgr = image.at<cv::Vec3b>(ii, jj);

            bgr[0] = finalz;
            bgr[1] = finalz;
            bgr[2] = finalz;
          }
        }
      }

      // Save results.
      std::cout << "Resolution: " << rows << "x" << cols << std::endl;
      std::cout << "Rows: " << diff1 << " Cols: " << diff2 << std::endl;
      std::cout << "Val min: " << m_val_min << " Val max: " << m_val_max << std::endl;
      std::cout << "Num cells: " << m_grid.size() << " : " << numcells << std::endl;

      imwrite("/Users/monet/Documents/lod/logs/buildings/image.jpg", image);
    }

    void apply_lsd() {

      Mat image = imread("/Users/monet/Documents/lod/logs/buildings/image.jpg", cv::IMREAD_GRAYSCALE);
      imwrite("/Users/monet/Documents/lod/logs/buildings/image-gscale.jpg", image);

      cv::Ptr<cv::LineSegmentDetector> ls = createLineSegmentDetector(
        cv::LSD_REFINE_STD, m_lsd_scale);

      double start = double(cv::getTickCount());
      vector<cv::Vec4f> lines_std;

      // Detect the lines.
      ls->detect(image, lines_std);
      double duration_ms = (double(cv::getTickCount()) - start) * 1000 / cv::getTickFrequency();
      std::cout << "It took " << duration_ms << " ms." << std::endl;

      // Show found lines.
      Mat drawn_lines(image);
      ls->drawSegments(drawn_lines, lines_std);
      imwrite("/Users/monet/Documents/lod/logs/buildings/result-lsd.jpg", drawn_lines);
    }

    void save_cluster(const std::string name) {
      
      std::vector<Point_3> points;
      points.reserve(m_cluster.size());
      for (const auto& item: m_cluster)
        points.push_back(item.input_point);
      const Color color(0, 0, 0);
      m_saver.export_points(points, color, name);
    }
  };

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_INTERNAL_CLOUD_TO_IMAGE_CONVERTER_H
