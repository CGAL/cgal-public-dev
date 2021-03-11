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

#ifndef CGAL_LEVELS_OF_DETAIL_INTERNAL_PLANAR_IMAGE_REGION_H
#define CGAL_LEVELS_OF_DETAIL_INTERNAL_PLANAR_IMAGE_REGION_H

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
  class Planar_image_region {

  public:
    using Traits = GeomTraits;
    using Pixels = std::vector<Pixel>;
    using Indices = std::vector<std::size_t>;

    Planar_image_region(
      const Pixels& pixels) :
    m_pixels(pixels)
    { }

    bool is_already_visited(
      const std::size_t, const std::size_t, const bool) const {
      return false;
    }

    bool is_part_of_region(
      const std::size_t from_index,
      const std::size_t to_index,
      const Indices& region) const {

      CGAL_precondition(region.size() > 0);
      const Pixel& from = m_pixels[from_index];
      const Pixel& to   = m_pixels[to_index];

      return (
        from.label == to.label && to.label != std::size_t(-1));
    }

    bool is_valid_region(const Indices& region) const {
      return region.size() >= 1;
    }

    void update(const Indices&) {
      // skipped!
    }

  private:
    const Pixels& m_pixels;
  };

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_INTERNAL_PLANAR_IMAGE_REGION_H
