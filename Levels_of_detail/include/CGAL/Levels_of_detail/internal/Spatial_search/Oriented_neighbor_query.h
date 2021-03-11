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

#ifndef CGAL_LEVELS_OF_DETAIL_INTERNAL_ORIENTED_NEIGHBOR_QUERY_H
#define CGAL_LEVELS_OF_DETAIL_INTERNAL_ORIENTED_NEIGHBOR_QUERY_H

// STL includes.
#include <vector>

// CGAL includes.
#include <CGAL/assertions.h>

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

  template<
  typename GeomTraits,
  typename Segment_wrapper>
  class Oriented_neighbor_query {

  public:

    using Traits = GeomTraits;
    using Indices = std::vector<std::size_t>;

    Oriented_neighbor_query(
      const std::vector<Segment_wrapper>& wrappers,
      const Indices& idx_map) :
    m_wrappers(wrappers),
    m_idx_map(idx_map)
    { }

    void operator()(
      const std::size_t query_index,
      Indices& neighbors) const {

      CGAL_assertion(query_index >= 0);
      CGAL_assertion(query_index < m_idx_map.size());

      neighbors.clear();
      const std::size_t real_index = m_idx_map[query_index];
      const auto& wrapper = m_wrappers[real_index];
      neighbors = wrapper.neighbors;
    }

  private:
    const std::vector<Segment_wrapper>& m_wrappers;
    const Indices& m_idx_map;
  };

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_INTERNAL_ORIENTED_NEIGHBOR_QUERY_H
