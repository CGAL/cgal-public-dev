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

#ifndef CGAL_URBAN_AREA_PROCESSING_UTILS_H
#define CGAL_URBAN_AREA_PROCESSING_UTILS_H

// #include <CGAL/license/Urban_area_processing.h>

// Internal includes.
#include <CGAL/Urban_area_processing/internal/Tools/Generic_point_extractor.h>
#include <CGAL/Urban_area_processing/internal/Contouring/Boundary_from_triangulation_2.h>
#include <CGAL/Urban_area_processing/internal/Contouring/Shortest_path_contouring_2.h>

// TODO:
// 1. Use input traits, now they are passed but never used.

namespace CGAL {
namespace Urban_area_processing {

  void triangle_mesh_to_point_cloud() {

  }

  template<
  typename GeomTraits,
  typename InputRange,
  typename OutputCondition,
  typename PointMap,
  typename OutputIterator>
  void filter_points(
    GeomTraits& traits,
    const InputRange& input_range,
    const OutputCondition& condition,
    const PointMap point_map,
    OutputIterator points) {

    using Point_extractor = internal::Generic_point_extractor<
    GeomTraits, InputRange, OutputCondition, PointMap>;

    const Point_extractor extractor(
      input_range, condition, point_map);
    extractor.extract(points);
  }

  void split_into_clusters() {

  }

  template<
  typename GeomTraits,
  typename InputRange,
  typename SegmentMap,
  typename OutputIterator>
  void merge_and_orient_segments(
    GeomTraits& traits,
    const InputRange& input_range,
    const SegmentMap segment_map,
    OutputIterator contours,
    const typename GeomTraits::FT scale,
    const typename GeomTraits::FT min_length_2) {

    using Shortest_path_contouring_2 = 
      internal::Shortest_path_contouring_2<GeomTraits, InputRange, SegmentMap>;
    
    Shortest_path_contouring_2 shortest(
      input_range, segment_map, scale, min_length_2, true);
    shortest.merge(contours);
  }

  void merge_and_orient_contours() {

  }

  void label_triangulation() {
    
  }

  template<
  typename GeomTraits,
  typename InputTriangulation,
  typename OutputIterator>
  void extract_boundary_with_holes(
    GeomTraits& traits,
    InputTriangulation& triangulation,
    OutputIterator boundaries) {

    using Boundary_extractor = 
      internal::Boundary_from_triangulation_2<GeomTraits, InputTriangulation>;

    Boundary_extractor extractor(triangulation);
    extractor.extract(boundaries);
  }

  void refine_triangulation() {

  }

} // Urban_area_processing
} // CGAL

#endif // CGAL_URBAN_AREA_PROCESSING_UTILS_H
