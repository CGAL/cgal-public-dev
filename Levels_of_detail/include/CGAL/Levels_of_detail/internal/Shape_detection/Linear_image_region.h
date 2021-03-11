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

#ifndef CGAL_LEVELS_OF_DETAIL_INTERNAL_LINEAR_IMAGE_REGION_H
#define CGAL_LEVELS_OF_DETAIL_INTERNAL_LINEAR_IMAGE_REGION_H

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
  typename Pixel>
  class Linear_image_region {

  public:
    using Traits = GeomTraits;
    using Indices = std::vector<std::size_t>;

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
      return region.size() >= 2;
    }

    void update(const Indices&) {
      // skipped!
    }
  };

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_INTERNAL_LINEAR_IMAGE_REGION_H
