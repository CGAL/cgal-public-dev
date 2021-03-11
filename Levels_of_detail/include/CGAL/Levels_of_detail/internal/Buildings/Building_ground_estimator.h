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

#ifndef CGAL_LEVELS_OF_DETAIL_INTERNAL_BUILDING_GROUND_ESTIMATOR_H
#define CGAL_LEVELS_OF_DETAIL_INTERNAL_BUILDING_GROUND_ESTIMATOR_H

// STL includes.
#include <vector>

// CGAL includes.
#include <CGAL/number_utils.h>

// Internal includes.
#include <CGAL/Levels_of_detail/internal/utils.h>
#include <CGAL/Levels_of_detail/internal/struct.h>

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

  template<
  typename GeomTraits,
  typename InputTriangulation>
  class Building_ground_estimator {

  public:
    using Traits = GeomTraits;
    using Triangulation = InputTriangulation;

    using FT = typename Traits::FT;
    using Point_3 = typename Traits::Point_3;

    using Approximate_face = internal::Partition_edge_3<Traits>;

    Building_ground_estimator(
      const Triangulation& base,
      const FT ground_z) :
    m_base(base),
    m_ground_z(ground_z)
    { }

    void estimate(
      Approximate_face& ground) const {

      FT minx = internal::max_value<FT>(), miny = internal::max_value<FT>();
      FT maxx = -internal::max_value<FT>(), maxy = -internal::max_value<FT>();

      for (auto fh = m_base.finite_faces_begin();
      fh != m_base.finite_faces_end(); ++fh) {
        for (std::size_t k = 0; k < 3; ++k) {
          const auto& p = fh->vertex(k)->point();
          minx = CGAL::min(minx, p.x()); miny = CGAL::min(miny, p.y());
          maxx = CGAL::max(maxx, p.x()); maxy = CGAL::max(maxy, p.y());
        }
      }

      ground.polygon.clear(); ground.polygon.resize(4);
      ground.polygon[0] = Point_3(minx, miny, m_ground_z);
      ground.polygon[1] = Point_3(maxx, miny, m_ground_z);
      ground.polygon[2] = Point_3(maxx, maxy, m_ground_z);
      ground.polygon[3] = Point_3(minx, maxy, m_ground_z);
    }

  private:
    const Triangulation& m_base;
    const FT m_ground_z;
  };

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_INTERNAL_BUILDING_GROUND_ESTIMATOR_H
