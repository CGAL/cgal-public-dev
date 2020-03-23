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
#include <CGAL/Urban_area_processing/internal/Generic_point_extractor.h>
#include <CGAL/Urban_area_processing/internal/Boundary_from_triangulation_2.h>

// TODO:
// 1. Use input traits, now they are passed but never used.

namespace CGAL {
namespace Urban_area_processing {

  void triangle_mesh_to_point_cloud() {

  }

  void export_urban_object_in_city_GML() {

  }

  template<
  typename GeomTraits,
  typename InputRange,
  typename OutputCondition,
  typename PointMap,
  typename OutputIterator>
  void extract_points(
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

  void merge_and_orient_segments() {

  }

  void merge_and_orient_contours() {

  }

  template<
  typename GeomTraits,
  typename InputTriangulation,
  typename OutputIterator>
  void extract_boundary_with_holes_from_triangulation(
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

  void extrude_triangulation() {

  }

  void mark_inside_outside_in_triangulation() {
    
  }

} // Urban_area_processing
} // CGAL

#endif // CGAL_URBAN_AREA_PROCESSING_UTILS_H
