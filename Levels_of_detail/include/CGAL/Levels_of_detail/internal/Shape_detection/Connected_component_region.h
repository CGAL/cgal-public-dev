// Copyright (c) 2019 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
//
//
// Author(s)     : Dmitry Anisimov, Simon Giraudot, Pierre Alliez, Florent Lafarge, and Andreas Fabri
//

#ifndef CGAL_LEVELS_OF_DETAIL_INTERNAL_SHAPE_DETECTION_CONNECTED_COMPONENT_REGION_H
#define CGAL_LEVELS_OF_DETAIL_INTERNAL_SHAPE_DETECTION_CONNECTED_COMPONENT_REGION_H

#include <CGAL/license/Levels_of_detail.h>

// STL includes.
#include <vector>

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

  class Connected_component_region {
  public:
    Connected_component_region(
      const std::size_t min_region_size) :
    m_min_region_size(min_region_size)
    { }

    bool is_already_visited(
      const std::size_t,
      const std::size_t query_index,
      const bool is_visited) const { return false; }

    bool is_part_of_region(
      const std::size_t,
      const std::size_t,
      const std::vector<std::size_t>&) const {
      return true;
    }

    bool is_valid_region(const std::vector<std::size_t>& region) const {
      return ( region.size() >= m_min_region_size );
    }

    void update(const std::vector<std::size_t>&) { }
  private:
    const std::size_t m_min_region_size;
  };

} // namespace internal
} // namespace Levels_of_detail
} // namespace CGAL

#endif // CGAL_LEVELS_OF_DETAIL_INTERNAL_SHAPE_DETECTION_CONNECTED_COMPONENT_REGION_H
