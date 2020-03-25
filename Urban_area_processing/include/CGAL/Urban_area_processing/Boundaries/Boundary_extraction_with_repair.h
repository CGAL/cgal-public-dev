// Copyright (c) 2020 SARL GeometryFactory (France).
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
// Author(s)     : Dmitry Anisimov, Simon Giraudot, Pierre Alliez, Florent Lafarge, and Andreas Fabri

#ifndef CGAL_URBAN_AREA_PROCESSING_BOUNDARY_EXTRACTION_WITH_REPAIR_H
#define CGAL_URBAN_AREA_PROCESSING_BOUNDARY_EXTRACTION_WITH_REPAIR_H

// #include <CGAL/license/Urban_area_processing.h>

// Boost includes.
#include <boost/function_output_iterator.hpp>

// CGAL includes.
#include <CGAL/assertions.h>
#include <CGAL/property_map.h>
#include <CGAL/Shape_detection/Region_growing/Region_growing.h>

// Internal includes.
#include <CGAL/Urban_area_processing/struct.h>
#include <CGAL/Urban_area_processing/internal/utils.h>
#include <CGAL/Urban_area_processing/internal/property_map.h>
#include <CGAL/Urban_area_processing/internal/Tools/Generic_point_extractor.h>
#include <CGAL/Urban_area_processing/internal/Tools/Extract_vertical_points_3.h>
#include <CGAL/Urban_area_processing/internal/Tools/Point_3_to_point_2_inserter.h>
#include <CGAL/Urban_area_processing/internal/Shape_detection/Estimate_normals_2.h>
#include <CGAL/Urban_area_processing/internal/Shape_detection/Estimate_normals_3.h>
#include <CGAL/Urban_area_processing/internal/Shape_detection/Sphere_neighbor_query.h>
#include <CGAL/Urban_area_processing/internal/Shape_detection/Least_squares_line_fit_region_2.h>
#include <CGAL/Urban_area_processing/internal/Shape_detection/Least_squares_line_fit_sorting_2.h>
#include <CGAL/Urban_area_processing/internal/Contouring/Shortest_path_contouring_2.h>
#include <CGAL/Urban_area_processing/internal/Contouring/Boundary_from_triangulation_2.h>

// Utils.
#include "../../../../test/Urban_area_processing/include/Saver.h"

namespace CGAL {
namespace Urban_area_processing {

  /*!
    \ingroup PkgUrbanAreaProcessingRefBoundaries

    \brief extracts an approximate closed contour, possibly with holes, 
    from the noisy point cloud with some missing boundary points.

    This class 
    - extracts all boundary points;
    - projects them onto the xy plane;
    - detects linear regions in this point set;
    - fit segments to the detected regions;
    - merges and orient segments into open contours;
    - merges all open contours into a closed outer boundary contour;
    - detects holes with respect to the outer contour.

    This class assumes that input point cloud is of varying density and possibly 
    with missing boundary points. But missing boundary regions should not exceed 
    more than 30% of the whole boundary.

    \tparam GeomTraits 
    must be a model of `Kernel`.

    \tparam InputRange
    must be a model of `ConstRange` whose iterator type is `RandomAccessIterator`.

    \tparam PointMap 
    must be an `LvaluePropertyMap` whose key type is the value type of the input 
    range and value type is `GeomTraits::Point_3`.
  */
  template<
  typename GeomTraits,
  typename InputRange,
  typename PointMap>
  class Boundary_extraction_with_repair {

  public:
    
    /// \cond SKIP_IN_MANUAL
    using Traits = GeomTraits;
    using Input_range = InputRange;
    using Point_map = PointMap;

    using FT = typename Traits::FT;
    using Point_2 = typename Traits::Point_2;
    using Vector_2 = typename Traits::Vector_2;
    using Vector_3 = typename Traits::Vector_3;
    using Segment_2 = typename Traits::Segment_2;
    using Line_2 = typename Traits::Line_2;

    using Indices = std::vector<std::size_t>;
    using Points_2 = std::vector<Point_2>;
    using Pair_item_2 = std::pair<Point_2, Vector_2>;
    using Pair_range_2 = std::vector<Pair_item_2>;
    using First_of_pair_map = CGAL::First_of_pair_property_map<Pair_item_2>;
    using Second_of_pair_map = CGAL::Second_of_pair_property_map<Pair_item_2>;
    using Segments_2 = std::vector<Segment_2>;

    using Sphere_neighbor_query_3 = internal::Sphere_neighbor_query<
      Traits, Input_range, Point_map>;
    using Estimate_normals_3 = internal::Estimate_normals_3<
      Traits, Input_range, Sphere_neighbor_query_3>;

    using Point_inserter = 
      internal::Point_3_to_point_2_inserter<Traits>;
    using Vertical_condition = 
      internal::Extract_vertical_points_3<Traits>;
    using Point_extractor = internal::Generic_point_extractor<
      Traits, Input_range, Vertical_condition, Point_map>;

    using Identity_map_2 = 
      CGAL::Identity_property_map<Point_2>;
    using Sphere_neighbor_query_2 = internal::Sphere_neighbor_query<
      Traits, Points_2, Identity_map_2>;

    using Estimate_normals_2 = internal::Estimate_normals_2<
      Traits, Points_2, Sphere_neighbor_query_2>;
    using LSLF_region = internal::Least_squares_line_fit_region_2<
      Traits, Pair_range_2, First_of_pair_map, Second_of_pair_map>;
    using LSLF_sorting = internal::Least_squares_line_fit_sorting_2<
      Traits, Points_2, Sphere_neighbor_query_2, Identity_map_2>;
    using Region_growing_2 = CGAL::Shape_detection::Region_growing<
      Points_2, Sphere_neighbor_query_2, LSLF_region, typename LSLF_sorting::Seed_map>;

    using Boundary_point_map_2 = 
      internal::Item_property_map<Points_2, Identity_map_2>;
    
    using Segment_map_2 = 
      CGAL::Identity_property_map<Segment_2>;
    using Shortest_path_contouring_2 = 
      internal::Shortest_path_contouring_2<Traits, Segments_2, Segment_map_2>;

    using Triangulation = Triangulation<Traits>;
    using Delaunay = typename Triangulation::Delaunay;
    using Boundary_extractor = 
      internal::Boundary_from_triangulation_2<Traits, Delaunay>;
    /// \endcond

    /// \name Initialization
    /// @{

    /*!
      \brief initializes all internal data structures.
      
      \param input_range
      an instance of `InputRange` with 3D points

      \param point_map
      an instance of `PointMap` that maps an item from `input_range` 
      to `GeomTraits::Point_3`

      \param scale
      a user-defined scale in meters

      \param noise
      a user-defined noise level in meters

      \param min_length_2
      a min length of the boundary segment in meters

      \param max_angle_2
      will be redefined later (in degrees)

      \param max_angle_3
      will be redefined later (in degrees)

      \pre `scale > 0`
      \pre `noise > 0`
      \pre `min_length_2 > 0`
      \pre `max_angle_2 > 0 && max_angle_2 < 90`
      \pre `max_angle_3 > 0 && max_angle_3 < 90`
    */
    Boundary_extraction_with_repair(
      const InputRange& input_range,
      const PointMap point_map,
      const FT scale,
      const FT noise,
      const FT min_length_2,
      const FT max_angle_2,
      const FT max_angle_3) : 
    m_input_range(input_range),
    m_point_map(point_map),
    m_scale(scale),
    m_noise(noise),
    m_min_length_2(min_length_2),
    m_max_angle_2(max_angle_2),
    m_max_angle_3(max_angle_3) { 
      
      CGAL_precondition(input_range.size() > 0);
    }

    /// @}

    /// \name Extraction
    /// @{

    /*!
      \brief extracts a set of boundary contours.

      \tparam OutputIterator 
      must be an output iterator whose value type is `std::vector< std::pair<Point_3, std::size_t> >`,
      where the first item in the pair is a point and second item is the contour index. 
      If the latter is `std::size_t(-1)` then this contour is outer, otherwise it is a hole
      and the stored index is the index of the corresponding outer contour.

      \param boundaries
      an output iterator with boundary contours
    */
    template<typename OutputIterator>
    void extract(OutputIterator boundaries) {

      std::cout << "* extracting boundary with repair... " << std::endl;
      std::vector<Point_2> boundary_points_2;
      extract_boundary_points_2(boundary_points_2);
      std::cout << "- boundary points are extracted: " << 
        boundary_points_2.size() << std::endl;

      save_boundary_points_2(boundary_points_2);
      std::cout << "- boundary points are saved" << std::endl;

      std::vector<Indices> wall_regions_2;
      extract_wall_regions_2(boundary_points_2, wall_regions_2);
      std::cout << "- wall regions are extracted: " << 
        wall_regions_2.size() << std::endl;
      
      save_wall_regions_2(boundary_points_2, wall_regions_2);
      std::cout << "- wall regions are saved" << std::endl;

      std::vector<Segment_2> segments;
      create_wall_segments(
        boundary_points_2, wall_regions_2, segments);
      std::cout << "- wall segments are created: " << 
        segments.size() << std::endl;

      save_wall_segments(segments, "wall_segments_original");
      std::cout << "- wall segments are saved " << std::endl;

      const std::size_t num_subcontours = merge_segments(segments);
      std::cout << "- wall subcontours are created: " << 
        num_subcontours << std::endl;
      
      save_wall_segments(segments, "wall_segments_merged");
      std::cout << "- wall subcontours are saved " << std::endl;

      Triangulation triangulation;
      close_outer_boundary(segments, triangulation);
      std::cout << "- outer boundary is closed " << std::endl;

      const std::size_t num_holes = add_holes(triangulation);
      std::cout << "- holes are added: " << num_holes << std::endl;

      /*
      Boundary_extractor extractor(triangulation.delaunay, false);
      extractor.extract(boundaries); 
      std::cout << "- boundaries and holes are extracted " << std::endl;
      */

      std::cout << std::endl;
    }

    /// @}

  private:
    const Input_range& m_input_range;
    const Point_map m_point_map;
    const FT m_scale;
    const FT m_noise;
    const FT m_min_length_2;
    const FT m_max_angle_2;
    const FT m_max_angle_3;

    void extract_boundary_points_2(
      std::vector<Point_2>& boundary_points_2) const {

      std::vector<Vector_3> normals;
      estimate_normals(normals);

      boundary_points_2.clear();
      const Vertical_condition vertical_condition(
        normals, m_max_angle_3);
      const Point_extractor extractor(
        m_input_range, vertical_condition, m_point_map);

      Point_inserter inserter(boundary_points_2);
      extractor.extract(
        boost::make_function_output_iterator(inserter));
    }

    void estimate_normals(
      std::vector<Vector_3>& normals) const {

      Sphere_neighbor_query_3 neighbor_query(
        m_input_range, m_scale, m_point_map);
      Estimate_normals_3 estimator(
        m_input_range, neighbor_query);
      estimator.get_normals(normals);
      CGAL_assertion(normals.size() == m_input_range.size());
    }

    void save_boundary_points_2(
      const std::vector<Point_2>& boundary_points_2) const {

      Saver<Traits> saver;
      saver.export_points(boundary_points_2,  
      "/Users/monet/Documents/gf/urban-area-processing/logs/boundary_points_2");
    }

    void extract_wall_regions_2(
      const std::vector<Point_2>& boundary_points_2,
      std::vector<Indices>& wall_regions_2) const {

      wall_regions_2.clear();
      Identity_map_2 identity_map_2;
      Sphere_neighbor_query_2 neighbor_query(
        boundary_points_2, m_scale, identity_map_2);

      std::vector<Vector_2> normals;
      Estimate_normals_2 estimator(
        boundary_points_2, neighbor_query);
      estimator.get_normals(normals);

      CGAL_assertion(boundary_points_2.size() == normals.size());
      Pair_range_2 range;
      range.reserve(boundary_points_2.size());
      for (std::size_t i = 0; i < boundary_points_2.size(); ++i)
        range.push_back(std::make_pair(boundary_points_2[i], normals[i]));

      First_of_pair_map point_map;
      Second_of_pair_map normal_map;
      LSLF_region region(
        range, 
        m_noise,
        m_max_angle_2,
        m_min_length_2,
        point_map,
        normal_map);

      LSLF_sorting sorting(
        boundary_points_2, neighbor_query, identity_map_2);
      sorting.sort();

      Region_growing_2 region_growing(
        boundary_points_2,
        neighbor_query,
        region,
        sorting.seed_map());
      region_growing.detect(std::back_inserter(wall_regions_2));
    }

    void save_wall_regions_2(
      const std::vector<Point_2>& boundary_points_2,
      const std::vector<Indices>& wall_regions_2) const {

      Saver<Traits> saver;
      saver.export_points(boundary_points_2, wall_regions_2,
      "/Users/monet/Documents/gf/urban-area-processing/logs/wall_regions_2");
    }

    void create_wall_segments(
      const std::vector<Point_2>& boundary_points_2,
      const std::vector<Indices>& wall_regions_2,
      std::vector<Segment_2>& segments) const {

      segments.clear();
      segments.reserve(wall_regions_2.size());

      Identity_map_2 identity_map_2;
      Boundary_point_map_2 point_map_2(
        boundary_points_2, identity_map_2);

      Indices indices;
      Line_2 line; Point_2 p, q; 
      for (const auto& wall_region : wall_regions_2) {
        indices.clear();
        for (std::size_t i = 0; i < wall_region.size(); ++i)
          indices.push_back(i);
        internal::line_from_points_2(
          wall_region, point_map_2, line);
        internal::boundary_points_on_line_2(
          wall_region, point_map_2, indices, line, p, q);
        segments.push_back(Segment_2(p, q));
      }
      CGAL_assertion(
        segments.size() == wall_regions_2.size());
    }

    void save_wall_segments(
      const std::vector<Segment_2>& segments,
      const std::string name) const {

      Saver<Traits> saver;
      saver.export_polylines(segments,
      "/Users/monet/Documents/gf/urban-area-processing/logs/" + name);
    }

    std::size_t merge_segments(
      std::vector<Segment_2>& segments) const {

      Segment_map_2 segment_map;
      std::vector<Segments_2> subcontours;
      Shortest_path_contouring_2 shortest(
        segments, segment_map, m_scale, m_min_length_2, false);
      shortest.merge(std::back_inserter(subcontours));

      segments.clear();
      for (const auto& contour : subcontours)
        for (const auto& segment : contour)
          segments.push_back(segment);
      return subcontours.size();
    }

    void close_outer_boundary(
      const std::vector<Segment_2>& segments,
      Triangulation& triangulation) const {

      
    }

    std::size_t add_holes(Triangulation& triangulation) const {

      return 0;
    }
  };

} // Urban_area_processing
} // CGAL

#endif // CGAL_URBAN_AREA_PROCESSING_BOUNDARY_EXTRACTION_WITH_REPAIR_H
