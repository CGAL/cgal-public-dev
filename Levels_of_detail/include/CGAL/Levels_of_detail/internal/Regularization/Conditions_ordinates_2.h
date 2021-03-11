// Copyright (c) 2019 GeometryFactory Sarl (France).
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
// Author(s)     : Jean-Philippe Bauchet, Florent Lafarge, Gennadii Sytov, Dmitry Anisimov
//

#ifndef CGAL_LEVELS_OF_DETAIL_INTERNAL_CONDITIONS_ORDINATES_2
#define CGAL_LEVELS_OF_DETAIL_INTERNAL_CONDITIONS_ORDINATES_2

#include <CGAL/license/Levels_of_detail.h>

#include <CGAL/Levels_of_detail/internal/utils.h>
#include <CGAL/Levels_of_detail/internal/Regularization/Segment_data_2.h>

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

  template<
    typename GeomTraits>
  class Conditions_ordinates_2 {

  public:
    using Traits = GeomTraits;
    using FT = typename GeomTraits::FT;
    using Segment_data = typename internal::Segment_data_2<Traits>;

    Conditions_ordinates_2() :
    m_moe(FT(1)) { }

    FT reference(const Segment_data & seg_data, const FT suffix) const {
    FT val = seg_data.m_reference_coordinates.y() + suffix;
    return val;
    }

    int group_index(const FT val, const FT val_j, const int g_index) const {
      int g_j = -1;
      if (CGAL::abs(val_j - val) < m_moe)
        g_j = g_index;
      return g_j;
    }

    FT get_margin_of_error() const {
      return m_moe;
    }

    void set_margin_of_error(const FT max_bound) {
      CGAL_precondition(max_bound > 0);
      m_moe = max_bound / FT(100);
    }

  private:
    FT m_moe;
  };

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_INTERNAL_CONDITIONS_ORDINATES_2
