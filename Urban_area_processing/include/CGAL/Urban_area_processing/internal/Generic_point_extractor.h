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

#ifndef CGAL_URBAN_AREA_PROCESSING_INTERNAL_GENERIC_POINT_EXTRACTOR_H
#define CGAL_URBAN_AREA_PROCESSING_INTERNAL_GENERIC_POINT_EXTRACTOR_H

// #include <CGAL/license/Urban_area_processing.h>

// STL includes.
#include <vector>
#include <utility>
#include <algorithm>

// CGAL includes.
#include <CGAL/assertions.h>

// TODO:
// Remove this class and substitute it by the corresponding boost filter.

namespace CGAL {
namespace Urban_area_processing {
namespace internal {

  template<
  typename GeomTraits,
  typename InputRange,
  typename OutputCondition,
  typename PointMap>
  class Generic_point_extractor {

  public:
    using Traits = GeomTraits;
    using Input_range = InputRange;
    using Output_condition = OutputCondition;
    using Point_map = PointMap;

    Generic_point_extractor(
      const Input_range& input_range,
      const Output_condition& condition,
      const Point_map point_map) :
    m_input_range(input_range),
    m_condition(condition),
    m_point_map(point_map) { 
      
      CGAL_precondition(input_range.size() > 0);
    }

    template<typename OutputIterator>
    void extract(OutputIterator points) const {
      for (std::size_t i = 0; i < m_input_range.size(); ++i) {
        const bool satisfied = m_condition.is_valid(i);
        if (satisfied)
          *(++points) = get(m_point_map, *(m_input_range.begin() + i));
      }
    }

  private:
    const Input_range& m_input_range;
    const Output_condition& m_condition;
    const Point_map m_point_map;
  };

} // internal
} // Urban_area_processing
} // CGAL

#endif // CGAL_URBAN_AREA_PROCESSING_INTERNAL_GENERIC_POINT_EXTRACTOR_H
