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

#ifndef CGAL_URBAN_AREA_PROCESSING_INTERNAL_POINT_3_TO_POINT_2_INSERTER_H
#define CGAL_URBAN_AREA_PROCESSING_INTERNAL_POINT_3_TO_POINT_2_INSERTER_H

// #include <CGAL/license/Urban_area_processing.h>

// STL includes.
#include <vector>
#include <utility>
#include <algorithm>

// CGAL includes.
#include <CGAL/assertions.h>

namespace CGAL {
namespace Urban_area_processing {
namespace internal {

  template<typename GeomTraits>
  class Point_3_to_point_2_inserter {
    
  public:
    using Traits = GeomTraits;
    
    using Point_2 = typename Traits::Point_2;
    using Point_3 = typename Traits::Point_3;
      
    using argument_type = Point_3;
    using result_type = Point_2;

    Point_3_to_point_2_inserter(
      std::vector<Point_2>& result) : 
    m_result(result) { 
      m_result.clear();
    }

    void operator()(const argument_type& arg) {
      m_result.push_back(
        internal::point_2_from_point_3(arg));
    }

  private:
    std::vector<Point_2>& m_result;
  };

} // internal
} // Urban_area_processing
} // CGAL

#endif // CGAL_URBAN_AREA_PROCESSING_INTERNAL_POINT_3_TO_POINT_2_INSERTER_H
