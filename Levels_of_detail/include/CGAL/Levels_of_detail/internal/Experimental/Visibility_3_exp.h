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

#ifndef CGAL_LEVELS_OF_DETAIL_INTERNAL_VISIBILITY_3_EXP_H
#define CGAL_LEVELS_OF_DETAIL_INTERNAL_VISIBILITY_3_EXP_H

// STL includes.
#include <vector>
#include <utility>

// CGAL includes.
#include <CGAL/assertions.h>

// Internal includes.
#include <CGAL/Levels_of_detail/enum.h>
#include <CGAL/Levels_of_detail/internal/struct.h>
#include <CGAL/Levels_of_detail/internal/utils.h>

// Testing.
#include "../../../../../test/Levels_of_detail/include/Saver.h"

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

  template<
  typename GeomTraits,
  typename InputRange,
  typename PointMap>
  class Visibility_3_exp {
			
  public:
    using Traits = GeomTraits;
    using Input_range = InputRange;
    using Point_map = PointMap;

    using FT = typename Traits::FT;
    using Point_2 = typename Traits::Point_2;
    using Point_3 = typename Traits::Point_3;
    using Vector_3 = typename Traits::Vector_3;
    using Line_3 = typename Traits::Line_3;
    using Plane_3 = typename Traits::Plane_3;
    using Triangle_2 = typename Traits::Triangle_2;

    using Indices = std::vector<std::size_t>;
    using Partition_3 = internal::Partition_3<Traits>;
    using Stats = std::pair<FT, FT>;
    using Face = typename Partition_3::Face;
    using Building = internal::Building<Traits>;
    
    using Triangulation = internal::Triangulation<Traits>;
    using Vertex_handle = typename Triangulation::Delaunay::Vertex_handle;

    using Saver = Saver<Traits>;

    Visibility_3_exp(
      const Input_range& input_range,
      const Point_map& point_map,
      const Building& building,
      const std::vector<Indices>& roof_points_3) :
    m_input_range(input_range),
    m_point_map(point_map),
    m_building(building),
    m_roof_points_3(roof_points_3)
    { }

    void compute(Partition_3& partition) const {
      
    }
    
  private:
    const Input_range& m_input_range;
    const Point_map& m_point_map;
    const Building& m_building;
    const std::vector<Indices>& m_roof_points_3;
  };

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_INTERNAL_VISIBILITY_3_EXP_H
