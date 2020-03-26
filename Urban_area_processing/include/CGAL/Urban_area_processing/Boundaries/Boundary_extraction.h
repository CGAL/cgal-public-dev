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

#ifndef CGAL_URBAN_AREA_PROCESSING_BOUNDARY_EXTRACTION_H
#define CGAL_URBAN_AREA_PROCESSING_BOUNDARY_EXTRACTION_H

// #include <CGAL/license/Urban_area_processing.h>

namespace CGAL {
namespace Urban_area_processing {

  /*!
    \ingroup PkgUrbanAreaProcessingRefBoundaries

    \brief extracts an approximate closed contour, possibly with holes.

    This class identifies a type of the input point cloud and applies 
    the best corresponding contouring algorithm.

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
  class Boundary_extraction {

  public:
    /// \cond SKIP_IN_MANUAL
    using Traits = GeomTraits;
    using Input_range = InputRange;
    using Point_map = PointMap;
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
    */
    Boundary_extraction(
      const InputRange& input_range,
      const PointMap point_map) : 
    m_input_range(input_range),
    m_point_map(point_map) { 

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

    }

    /// @}

  private:
    const Input_range& m_input_range;
    const Point_map m_point_map;
  };

} // Urban_area_processing
} // CGAL

#endif // CGAL_URBAN_AREA_PROCESSING_BOUNDARY_EXTRACTION_H
