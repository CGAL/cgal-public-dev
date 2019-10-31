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

#ifndef CGAL_LEVELS_OF_DETAIL_INTERNAL_ROOF_VISIBILITY_3_H
#define CGAL_LEVELS_OF_DETAIL_INTERNAL_ROOF_VISIBILITY_3_H

// STL includes.
#include <vector>
#include <utility>
#include <iostream>
#include <fstream>
#include <memory>
#include <algorithm>

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

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

  template<typename GeomTraits>
  class Roof_visibility_3 {
			
  public:
    using Traits = GeomTraits;

    using FT = typename Traits::FT;
    using Point_2 = typename Traits::Point_2;
    using Point_3 = typename Traits::Point_3;
    using Triangle_2 = typename Traits::Triangle_2;

    using Indices = std::vector<std::size_t>;
    using Partition_2 = internal::Partition_2<Traits>;
    using Partition_3 = internal::Partition_3<Traits>;
    using Face_2 = typename Partition_2::Face;
    using Face_3 = typename Partition_3::Face;

    Roof_visibility_3(
      const Partition_2& partition_2) :
    m_partition_2(partition_2),
    m_num_samples(100),
    m_random(0)
    { }

    void compute(Partition_3& partition_3) {

    }
    
  private:
    const Partition_2& m_partition_2;
    const std::size_t m_num_samples;
    Random m_random;
  };

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_INTERNAL_ROOF_VISIBILITY_3_H
