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

#ifndef CGAL_LEVELS_OF_DETAIL_INTERNAL_VISIBILITY_2_H
#define CGAL_LEVELS_OF_DETAIL_INTERNAL_VISIBILITY_2_H

// STL includes.
#include <vector>
#include <utility>

// CGAL includes.
#include <CGAL/assertions.h>
#include <CGAL/property_map.h>

// LOD includes.
#include <CGAL/Levels_of_detail/enum.h>

// Internal includes.
#include <CGAL/Levels_of_detail/internal/utils.h>
#include <CGAL/Levels_of_detail/internal/struct.h>

// Testing.
#include "../../../../../test/Levels_of_detail/include/Saver.h"

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

template<
typename GeomTraits,
typename InputRange,
typename PointMap2>
class Visibility_2 {

public:
  using Traits = GeomTraits;
  using Input_range = InputRange;
  using Point_map_2 = PointMap2;

  using FT = typename Traits::FT;
  using Point_2 = typename Traits::Point_2;
  using Point_3 = typename Traits::Point_3;
  using Triangle_2 = typename Traits::Triangle_2;

  using Points_2 = std::vector<Point_2>;
  using Indices = std::vector<std::size_t>;

  using Triangulation = typename internal::Triangulation<Traits>::Delaunay;
  using Location_type = typename Triangulation::Locate_type;

  using Partition_2 = internal::Partition_2<Traits>;

  using Saver = Saver<Traits>;
  using Color = CGAL::Color;

  Visibility_2(
    const Input_range& input_range,
    const Point_map_2& point_map_2,
    const FT threshold = FT(1) / FT(2)) :
  m_input_range(input_range),
  m_point_map_2(point_map_2),
  m_threshold(threshold)
  { }

  void compute(Partition_2& partition_2) {
    save_input_range("/Users/monet/Documents/lod/logs/buildings/tmp/visibility_points.jpg");

    
  }

private:
  const Input_range& m_input_range;
  const Point_map_2& m_point_map_2;
  const FT m_threshold;

  Saver m_saver;

  void save_input_range(const std::string name) {

    std::vector<Point_3> points;
    points.reserve(m_input_range.size());
    for (const auto& item : m_input_range) {
      const auto& p = get(m_point_map_2, item);
      points.push_back(Point_3(p.x(), p.y(), FT(0)));
    }
      
    const Color color(0, 0, 0);
    m_saver.export_points(points, color, name);
  }
};

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_INTERNAL_VISIBILITY_2_H
