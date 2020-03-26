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

  /// \ingroup PkgUrbanAreaProcessingRef
  /*! 
    \brief converts a triangle mesh into a point cloud by applying a Monte Carlo 
    sampling to each of its faces.

    \tparam GeomTraits 
    must be a model of `Kernel`.

    \tparam SurfaceMesh 
    must be a model of `FaceListGraph`.

    \tparam OutputIterator 
    must be an output iterator whose value type is `GeomTraits::Point_3`.

    \param traits 
    an instance of `GeomTraits`

    \param mesh 
    an instance of `SurfaceMesh` with a triangle mesh

    \param point_cloud
    an output iterator with the computed point cloud
  */ 
  template<
  typename GeomTraits,
  typename SurfaceMesh,
  typename OutputIterator>
  void triangle_mesh_to_point_cloud(
    const GeomTraits& traits, 
    const SurfaceMesh& mesh,
    OutputIterator point_cloud) {

  }

  /// \ingroup PkgUrbanAreaProcessingRef
  /*! 
    \brief given one of the available filters, this function extracts points with a specific 
    property from the input point cloud.

    \tparam GeomTraits 
    must be a model of `Kernel`.

    \tparam InputRange 
    must be a model of `ConstRange` whose iterator type is `RandomAccessIterator`.

    \tparam OutputCondition
    a set of requirements that all returned points should satisfy.

    \tparam PointMap 
    must be an `LvaluePropertyMap` whose key type is the value type of the input 
    range and value type is `GeomTraits::Point_3`.

    \tparam OutputIterator 
    must be an output iterator whose value type is `GeomTraits::Point_3`.

    \param traits 
    an instance of `GeomTraits`

    \param input_range
    an instance of `InputRange` with 3D points

    \param condition
    an instance of `OutputCondition`

    \param point_map
    an instance of `PointMap` that maps an item from `input_range` 
    to `GeomTraits::Point_3`

    \param points
    an output iterator with the extracted points
  */ 
  template<
  typename GeomTraits,
  typename InputRange,
  typename OutputCondition,
  typename PointMap,
  typename OutputIterator>
  void filter_points(
    const GeomTraits& traits,
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

  /// \ingroup PkgUrbanAreaProcessingRef
  /*! 
    \brief splits a point cloud into separate clusters, for example multiple 
    buildings into separate buildings.

    \tparam GeomTraits 
    must be a model of `Kernel`.

    \tparam InputRange 
    must be a model of `ConstRange` whose iterator type is `RandomAccessIterator`.

    \tparam PointMap 
    must be an `LvaluePropertyMap` whose key type is the value type of the input 
    range and value type is `GeomTraits::Point_3`.

    \tparam OutputIterator 
    must be an output iterator whose value type is `std::vector<GeomTraits::Point_3>`.

    \param traits 
    an instance of `GeomTraits`

    \param input_range
    an instance of `InputRange` with 3D points

    \param point_map
    an instance of `PointMap` that maps an item from `input_range` 
    to `GeomTraits::Point_3`

    \param clusters
    an output iterator with the clusters
  */
  template<
  typename GeomTraits,
  typename InputRange,
  typename PointMap,
  typename OutputIterator>
  void split_into_clusters(
    const GeomTraits& traits,
    const InputRange& input_range,
    const PointMap point_map,
    OutputIterator clusters) {

  }

  /// \ingroup PkgUrbanAreaProcessingRef
  /*! 
    \brief given a set of disconnected segments, this function enables to connect 
    and orient them into a set of open contours.

    \tparam GeomTraits 
    must be a model of `Kernel`.

    \tparam InputRange 
    must be a model of `ConstRange` whose iterator type is `RandomAccessIterator`.

    \tparam SegmentMap 
    must be an `LvaluePropertyMap` whose key type is the value type of the input 
    range and value type is `GeomTraits::Segment_2`.

    \tparam OutputIterator 
    must be an output iterator whose value type is `std::vector<GeomTraits::Segment_2>`,

    \param traits 
    an instance of `GeomTraits`

    \param input_range
    an instance of `InputRange` with 2D segments

    \param segment_map
    an instance of `SegmentMap` that maps an item from `input_range` 
    to `GeomTraits::Segment_2`

    \param contours
    an output iterator with the found open contours

    \param scale 
    a user-defined scale

    \param min_length_2
    a user-defined min segment length
  */
  template<
  typename GeomTraits,
  typename InputRange,
  typename SegmentMap,
  typename OutputIterator>
  void merge_and_orient_segments(
    const GeomTraits& traits,
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

  /// \ingroup PkgUrbanAreaProcessingRef
  /*! 
    \brief given a set of open contours, this function tries to close and connect 
    them into a closed contour.

    \tparam GeomTraits 
    must be a model of `Kernel`.

    \tparam InputRange 
    must be a model of `ConstRange` whose iterator type is `RandomAccessIterator`.

    \tparam SegmentMap 
    must be an `LvaluePropertyMap` whose key type is the value type of the input 
    range and value type is `std::pair<GeomTraits::Segment_2, std::size_t>`, where 
    the first item in the pair is the segment and its second item is the index of the 
    contour this segment belongs to.

    \tparam OutputIterator 
    must be an output iterator whose value type is `std::vector<GeomTraits::Point_2>`,

    \param traits 
    an instance of `GeomTraits`

    \param input_range
    an instance of `InputRange` with 2D segments

    \param segment_map
    an instance of `SegmentMap` that maps an item from `input_range` 
    to `std::pair<GeomTraits::Segment_2, std::size_t>`

    \param contour
    an output iterator with the found closed contour
  */
  template<
  typename GeomTraits,
  typename InputRange,
  typename SegmentMap,
  typename OutputIterator>
  void merge_and_orient_contours(
    const GeomTraits& traits,
    const InputRange& input_range,
    const SegmentMap segment_map,
    OutputIterator contour) {

  }

  /// \ingroup PkgUrbanAreaProcessingRef
  /*! 
    \brief given a user-defined criteria, this function labels all triangulation 
    faces either inside or outside with respect to the urban object.

    \tparam GeomTraits 
    must be a model of `Kernel`.

    \tparam InputTriangulation 
    must be a model of `CGAL::Triangulation_2`.

    \param traits 
    an instance of `GeomTraits`

    \param triangulation
    an instance of `InputTriangulation` with 2D triangulation
  */
  template<
  typename GeomTraits,
  typename InputTriangulation>
  void label_triangulation(
    const GeomTraits& traits,
    InputTriangulation& triangulation) {
    
  }

  /// \ingroup PkgUrbanAreaProcessingRef
  /*! 
    \brief given a triangulation with face labels, that indicate either this face 
    belongs to the object interior or exterior, this function extracts a set of boundaries 
    from this triangulation.

    \tparam GeomTraits 
    must be a model of `Kernel`.

    \tparam InputTriangulation 
    must be a model of `CGAL::Triangulation_2`.

    \tparam OutputIterator 
    must be an output iterator whose value type is `std::vector< std::pair<Point_2, std::size_t> >`,
    where the first item in the pair is a point and second item is the contour index. 
    If the latter is `std::size_t(-1)` then this contour is outer, otherwise it is a hole
    and the stored index is the index of the corresponding outer contour.

    \param traits 
    an instance of `GeomTraits`

    \param triangulation
    an instance of `InputTriangulation` with 2D triangulation

    \param boundaries
    an output iterator with boundary contours
  */
  template<
  typename GeomTraits,
  typename InputTriangulation,
  typename OutputIterator>
  void extract_boundary_with_holes(
    const GeomTraits& traits,
    InputTriangulation& triangulation,
    OutputIterator boundaries) {

    using Boundary_extractor = 
      internal::Boundary_from_triangulation_2<GeomTraits, InputTriangulation>;

    Boundary_extractor extractor(triangulation);
    extractor.extract(boundaries);
  }

  /// \ingroup PkgUrbanAreaProcessingRef
  /*! 
    \brief refines triangulation from the point cloud with respect to the local 
    density of points.

    \tparam GeomTraits 
    must be a model of `Kernel`.

    \tparam InputTriangulation 
    must be a model of `CGAL::Triangulation_2`.

    \param traits 
    an instance of `GeomTraits`

    \param triangulation
    an instance of `InputTriangulation` with 2D triangulation
  */
  template<
  typename GeomTraits,
  typename InputTriangulation>
  void refine_triangulation(
    const GeomTraits& traits,
    InputTriangulation& triangulation) {

  }

} // Urban_area_processing
} // CGAL

#endif // CGAL_URBAN_AREA_PROCESSING_UTILS_H
