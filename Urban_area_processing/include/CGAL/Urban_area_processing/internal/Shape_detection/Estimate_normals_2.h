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

#ifndef CGAL_URBAN_AREA_PROCESSING_INTERNAL_ESTIMATE_NORMALS_2_H
#define CGAL_URBAN_AREA_PROCESSING_INTERNAL_ESTIMATE_NORMALS_2_H

// #include <CGAL/license/Urban_area_processing.h>

// STL includes.
#include <vector>

// CGAL includes.
#include <CGAL/assertions.h>

// TODO:
// Reimplement this class using CGAL normal 2 estimator options! Is that possible?

namespace CGAL {
namespace Urban_area_processing {
namespace internal {

  template<
  typename GeomTraits,
  typename InputRange,
  typename NeighborQuery>
  class Estimate_normals_2 {

  public:
    using Traits = GeomTraits;
    using Input_range = InputRange;
    using Neighbor_query = NeighborQuery;
    
    using FT = typename Traits::FT;
    using Vector_2 = typename Traits::Vector_2;
    using Line_2 = typename Traits::Line_2;

    using Indices = std::vector<std::size_t>;

    Estimate_normals_2(
      const Input_range& input_range,
      const Neighbor_query& neighbor_query) : 
    m_input_range(input_range),
    m_neighbor_query(neighbor_query) {
      
      CGAL_precondition(input_range.size() > 0);
    }

    void get_normals(
      std::vector<Vector_2>& normals) const {
                  
      normals.clear();
      normals.reserve(m_input_range.size());

      Indices neighbors;
      Line_2 line; Vector_2 normal;
      for (std::size_t i = 0; i < m_input_range.size(); ++i) {

        neighbors.clear();
        m_neighbor_query(i, neighbors);
        internal::line_from_points_2(
          neighbors, m_neighbor_query.point_map(), line);

        normal = line.to_vector();
        normal = normal.perpendicular(CGAL::COUNTERCLOCKWISE);

        const FT normal_length = internal::vector_length(normal);
        CGAL_assertion(normal_length > FT(0));
        normal /= normal_length;
        normals.push_back(normal);
      }
      CGAL_assertion(normals.size() == m_input_range.size());
    }

  private:
    const Input_range& m_input_range;
    const Neighbor_query& m_neighbor_query;
  };

} // namespace internal
} // namespace Urban_area_processing
} // namespace CGAL

#endif // CGAL_URBAN_AREA_PROCESSING_INTERNAL_ESTIMATE_NORMALS_2_H
