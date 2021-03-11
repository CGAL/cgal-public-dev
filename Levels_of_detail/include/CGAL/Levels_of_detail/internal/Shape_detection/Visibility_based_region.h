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

#ifndef CGAL_LEVELS_OF_DETAIL_INTERNAL_VISIBILITY_BASED_REGION_H
#define CGAL_LEVELS_OF_DETAIL_INTERNAL_VISIBILITY_BASED_REGION_H

// STL includes.
#include <vector>

// CGAL includes.
#include <CGAL/assertions.h>

// Internal includes.
#include <CGAL/Levels_of_detail/enum.h>

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

  template<
  typename GeomTraits,
  typename InputRange>
  class Visibility_based_region {

  public:
    using Traits = GeomTraits;
    using Input_range = InputRange;

    Visibility_based_region(
      const Input_range& input_range) :
    m_input_range(input_range)
    { }

    bool is_already_visited(
      const std::size_t,
      const std::size_t query_index,
      const bool is_visited) const { return false; }

    bool is_part_of_region(
      const std::size_t from_index,
      const std::size_t to_index,
      const std::vector<std::size_t>& region) const {

      if (region.size() == 1) {

        CGAL_precondition(region[0] >= 0);
        CGAL_precondition(region[0] < m_input_range.size());
        const auto& item = *(m_input_range.begin() + region[0]);
        if (item.visibility == Visibility_label::OUTSIDE)
          return false;
      }

      CGAL_precondition(
        from_index >= 0 &&
        to_index >= 0);
      CGAL_precondition(
        from_index < m_input_range.size() &&
        to_index < m_input_range.size());

      // Constrained-based.
      // const auto& from_item = *(m_input_range.begin() + from_index);
      // const auto& constr = from_item.constraints;
      // CGAL_assertion(constr.find(to_index) != constr.end());
      // if (constr.at(to_index))
        // return false;

      // Visibility-based.
      const auto& to_item = *(m_input_range.begin() + to_index);
      return to_item.visibility == Visibility_label::INSIDE;
    }

    bool is_valid_region(const std::vector<std::size_t>& region) const {

      if (region.size() == 1) {

        CGAL_precondition(region[0] >= 0);
        CGAL_precondition(region[0] < m_input_range.size());
        const auto& item = *(m_input_range.begin() + region[0]);
        if (item.visibility == Visibility_label::OUTSIDE)
          return false;
      }
      return region.size() >= 1;
    }

    void update(const std::vector<std::size_t>&) {
      // skipped!
    }

  private:
    const Input_range& m_input_range;
  };

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_INTERNAL_VISIBILITY_BASED_REGION_H
