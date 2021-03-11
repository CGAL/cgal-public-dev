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

#ifndef CGAL_LEVELS_OF_DETAIL_INTERNAL_SEGMENT_DATA_2
#define CGAL_LEVELS_OF_DETAIL_INTERNAL_SEGMENT_DATA_2

#include <CGAL/license/Levels_of_detail.h>

#include <CGAL/Levels_of_detail/internal/utils.h>

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

  template<
  typename GeomTraits>
  struct Segment_data_2 {

  public:
    using Traits = GeomTraits;
    using Segment = typename GeomTraits::Segment_2;
    using Point = typename GeomTraits::Point_2;
    using Vector  = typename GeomTraits::Vector_2;
    using FT = typename GeomTraits::FT;

    const Segment& m_segment;
    const std::size_t m_index;
    Vector  m_direction;
    FT      m_orientation;
    Point   m_barycentre;
    FT      m_length;
    Point   m_reference_coordinates;
    FT      m_a;
    FT      m_b;
    FT      m_c;

    Segment_data_2(
      const Segment& segment,
      const std::size_t index):
    m_segment(segment),
    m_index(index) {

      m_barycentre = middle_point_2(m_segment.source(), m_segment.target());
      m_length = static_cast<FT>(CGAL::sqrt(CGAL::to_double(m_segment.squared_length())));

      m_direction = internal::compute_direction(m_segment);
      m_orientation = internal::compute_orientation(m_direction);

      m_a = -static_cast<FT>(sin(CGAL::to_double(m_orientation * static_cast<FT>(CGAL_PI) / FT(180))));
      m_b =  static_cast<FT>(cos(CGAL::to_double(m_orientation * static_cast<FT>(CGAL_PI) / FT(180))));
      m_c = -m_a * m_barycentre.x() - m_b * m_barycentre.y();
    }
  };

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_INTERNAL_SEGMENT_DATA_2
