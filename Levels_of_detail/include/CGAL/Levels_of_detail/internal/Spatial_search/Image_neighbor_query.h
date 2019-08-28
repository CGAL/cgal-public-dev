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

#ifndef CGAL_LEVELS_OF_DETAIL_INTERNAL_IMAGE_NEIGHBOR_QUERY_H
#define CGAL_LEVELS_OF_DETAIL_INTERNAL_IMAGE_NEIGHBOR_QUERY_H

// STL includes.
#include <vector>

// CGAL includes.
#include <CGAL/assertions.h>

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

  template<
  typename GeomTraits,
  typename InputRange,
  typename ExtraRange,
  typename ExtraMap>
  class Image_neighbor_query {

  public:
    using Traits = GeomTraits;
    using Input_range = InputRange;
    using Extra_range = ExtraRange;
    using Extra_map = ExtraMap;

    using Index_to_point_map = 
    internal::Item_property_map<Extra_range, Extra_map>;

    Image_neighbor_query(
      const Input_range& input_range,
      const Extra_range& extra_range,
      const Extra_map& extra_map) :
    m_input_range(input_range),
    m_extra_range(extra_range),
    m_extra_map(extra_map),
    m_index_to_point_map(m_extra_range, m_extra_map)
    { }

    void operator()(
      const std::size_t query_index,
      std::vector<std::size_t>& neighbors) const {

      CGAL_assertion(query_index >= 0);
      CGAL_assertion(query_index < m_input_range.size());

      neighbors.clear();
      const auto& item = *(m_input_range.begin() + query_index);
      item.get_neighbors(neighbors);
    }

    const Index_to_point_map& point_map() const {
      return m_index_to_point_map;
    }

  private:
    const Input_range& m_input_range;
    const Extra_range& m_extra_range;
    const Extra_map& m_extra_map;

    const Index_to_point_map m_index_to_point_map;
  };

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_INTERNAL_IMAGE_NEIGHBOR_QUERY_H
