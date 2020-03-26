// All rights reserved.
// Copyright (c) 2020 SARL GeometryFactory (France).
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

#ifndef CGAL_URBAN_AREA_PROCESSING_INTERNAL_MERGE_CONTOURS_2_H
#define CGAL_URBAN_AREA_PROCESSING_INTERNAL_MERGE_CONTOURS_2_H

// #include <CGAL/license/Urban_area_processing.h>

// STL includes.
#include <vector>
#include <utility>
#include <algorithm>

// CGAL includes.
#include <CGAL/assertions.h>
#include <CGAL/property_map.h>

// Internal includes.
#include <CGAL/Urban_area_processing/internal/utils.h>

// TODO:
// Should I remove alpha shape computation at the beginning of the algorithm?

namespace CGAL {
namespace Urban_area_processing {
namespace internal {

  template<
  typename GeomTraits,
  typename InputRange,
  typename SegmentMap>
  class Merge_contours_2 {
  
  public:
    using Traits = GeomTraits;
    using Input_range = InputRange;
    using Segment_map = SegmentMap;

    using FT = typename Traits::FT;
    using Point_2 = typename Traits::Point_2;
    using Segment_2 = typename Traits::Segment_2;

    Merge_contours_2(
      const Input_range& input_range,
      const Segment_map segment_map,
      const FT scale,
      const FT noise,
      const bool verbose = true) :
    m_input_range(input_range),
    m_segment_map(segment_map),
    m_scale(scale),
    m_noise(noise),
    m_verbose(verbose) { 

      CGAL_precondition(input_range.size() > 0);
    }

    template<typename Triangulation>
    void merge(Triangulation& triangulation) const {

    }

  private:
    const Input_range& m_input_range;
    const Segment_map m_segment_map;
    const FT m_scale;
    const FT m_noise;
    const bool m_verbose;

  };

} // internal
} // Urban_area_processing
} // CGAL

#endif // CGAL_URBAN_AREA_PROCESSING_INTERNAL_MERGE_CONTOURS_2_H
