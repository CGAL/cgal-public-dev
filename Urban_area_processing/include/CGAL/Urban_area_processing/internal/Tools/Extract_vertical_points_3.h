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

#ifndef CGAL_URBAN_AREA_PROCESSING_INTERNAL_EXTRACT_VERTICAL_POINTS_3_H
#define CGAL_URBAN_AREA_PROCESSING_INTERNAL_EXTRACT_VERTICAL_POINTS_3_H

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
  class Extract_vertical_points_3 {

  public:
    using Traits = GeomTraits;

    using FT = typename Traits::FT;
    using Vector_3 = typename Traits::Vector_3;

    Extract_vertical_points_3(
      const std::vector<Vector_3>& normals,
      const FT max_angle_3) :
    m_normals(normals),
    m_max_angle_3(max_angle_3),
    m_vertical(
      Vector_3(FT(0), FT(0), FT(1))) {

      CGAL_precondition(normals.size() > 0);
      CGAL_precondition(max_angle_3 >= FT(0) && max_angle_3 <= FT(90));
    }

    bool is_valid(
      const std::size_t query_index) const {
      
      const auto& normal = m_normals[query_index];
      FT angle_3 = internal::angle_3d(normal, m_vertical);
      if (angle_3 > FT(90)) angle_3 = FT(180) - angle_3;
      angle_3 = FT(90) - angle_3;
      return angle_3 <= m_max_angle_3;
    }

  private:
    const std::vector<Vector_3>& m_normals;
    const FT m_max_angle_3;
    const Vector_3 m_vertical;
  };

} // internal
} // Urban_area_processing
} // CGAL

#endif // CGAL_URBAN_AREA_PROCESSING_INTERNAL_EXTRACT_VERTICAL_POINTS_3_H
