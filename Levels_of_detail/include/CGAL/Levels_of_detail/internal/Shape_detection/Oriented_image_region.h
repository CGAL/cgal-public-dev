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

#ifndef CGAL_LEVELS_OF_DETAIL_INTERNAL_ORIENTED_IMAGE_REGION_H
#define CGAL_LEVELS_OF_DETAIL_INTERNAL_ORIENTED_IMAGE_REGION_H

// STL includes.
#include <map>
#include <vector>
#include <utility>

// CGAL includes.
#include <CGAL/assertions.h>

// Internal includes.
#include <CGAL/Levels_of_detail/internal/utils.h>
#include <CGAL/Levels_of_detail/internal/struct.h>

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

  template<
  typename GeomTraits,
  typename Segment_wrapper>
  class Oriented_image_region {

  public:
    using Traits = GeomTraits;
    using Indices = std::vector<std::size_t>;

    Oriented_image_region(
      std::vector<Segment_wrapper>& wrappers,
      const Indices& idx_map,
      const std::size_t min_num = 2) :
    m_wrappers(wrappers),
    m_idx_map(idx_map),
    m_min_num(min_num)
    { }

    bool is_already_visited(
      const std::size_t, const std::size_t, const bool) const {
      return false;
    }

    bool is_part_of_region(
      const std::size_t, const std::size_t,
      const Indices& region) const {

      CGAL_precondition(region.size() > 0);
      return true;
    }

    bool is_valid_region(const Indices& region) const {
      return region.size() >= m_min_num;
    }

    void update(const Indices& region) {
      if (region.size() == 1) {
        const std::size_t real_index = m_idx_map[region[0]];
        auto& wrapper = m_wrappers[real_index];
        auto& neighbors = wrapper.neighbors;
        if (neighbors.size() >= 2) {
          for (std::size_t i = 0; i < neighbors.size(); ++i) {
            neighbors.pop_back();
            if (neighbors.size() == 1) break;
          }
        }
      }
    }

  private:
    std::vector<Segment_wrapper>& m_wrappers;
    const Indices& m_idx_map;
    const std::size_t m_min_num;
  };

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_INTERNAL_ORIENTED_IMAGE_REGION_H
