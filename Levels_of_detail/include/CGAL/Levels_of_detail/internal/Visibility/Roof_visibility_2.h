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

#ifndef CGAL_LEVELS_OF_DETAIL_INTERNAL_ROOF_VISIBILITY_2_H
#define CGAL_LEVELS_OF_DETAIL_INTERNAL_ROOF_VISIBILITY_2_H

// STL includes.
#include <vector>
#include <utility>
#include <iostream>
#include <fstream>
#include <memory>

// CGAL includes.
#include <CGAL/assertions.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/Random.h>

// Internal includes.
#include <CGAL/Levels_of_detail/enum.h>
#include <CGAL/Levels_of_detail/internal/struct.h>
#include <CGAL/Levels_of_detail/internal/utils.h>

// Spatial search.
#include <CGAL/Levels_of_detail/internal/Spatial_search/K_neighbor_query.h>

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

  template<typename GeomTraits>
  class Roof_visibility_2 {
			
  public:
    using Traits = GeomTraits;

    using FT = typename Traits::FT;
    using Point_2 = typename Traits::Point_2;
    using Point_3 = typename Traits::Point_3;
    using Vector_3 = typename Traits::Vector_3;
    using Plane_3 = typename Traits::Plane_3;

    using Indices = std::vector<std::size_t>;
    using Partition_2 = internal::Partition_2<Traits>;
    using Stats = std::pair<FT, FT>;
    using Face = typename Partition_2::Face;
    using Building = internal::Building<Traits>;

    using Generator = CGAL::Random_points_in_triangle_2<Point_2>;

    using Pair = std::pair<Point_2, FT>;
    using Point_map_2 = CGAL::First_of_pair_property_map<Pair>;
    using K_neighbor_query = internal::K_neighbor_query<Traits, std::vector<Pair>, Point_map_2>;
    using Random = CGAL::Random;

    Roof_visibility_2(
      const std::vector<Point_3>& input_range,
      const Building& building,
      const std::vector<Indices>& roof_points_3) :
    m_input_range(input_range),
    m_building(building),
    m_roof_points_3(roof_points_3),
    m_num_samples(100), // num samples per triangle
    m_k(6),
    m_random(0)
    { }

    void compute(Partition_2& partition) {
      
      const auto& ref = m_building.base1.triangulation.delaunay;

      for (auto& face : partition.faces) {
        const auto& tri = face.base.delaunay;
        const auto fh = tri.finite_faces_begin(); 

        const auto& p0 = fh->vertex(0)->point();
        const auto& p1 = fh->vertex(1)->point();
        const auto& p2 = fh->vertex(2)->point();

        const FT x = (p0.x() + p1.x() + p2.x()) / FT(3);
        const FT y = (p0.y() + p1.y() + p2.y()) / FT(3);

        const Point_2 b = Point_2(x, y);
        const auto bh = ref.locate(b);
        if (bh->info().tagged) {
          face.visibility = Visibility_label::INSIDE;
          face.inside = FT(1); face.outside = FT(0);
        } else {
          face.visibility = Visibility_label::OUTSIDE;
          face.inside = FT(0); face.outside = FT(1);
        }
      }
    }
    
  private:
    const std::vector<Point_3>& m_input_range;
    const Building& m_building;
    const std::vector<Indices>& m_roof_points_3;
    const std::size_t m_num_samples;
    std::vector<Point_3> m_samples;
    const std::size_t m_k;
    Random m_random;

    std::vector<std::size_t> m_roof_indices;
    std::vector< std::pair<Point_2, FT> > m_queries;
    Point_map_2 m_point_map_2;
    std::shared_ptr<K_neighbor_query> m_neighbor_query_ptr;
  };

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_INTERNAL_ROOF_VISIBILITY_2_H
