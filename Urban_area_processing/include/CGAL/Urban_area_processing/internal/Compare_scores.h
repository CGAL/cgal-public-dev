// All rights reserved.
// Copyright (c) 2020 SARL GeometryFactory (France).
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
// Author(s)     : Dmitry Anisimov, Simon Giraudot, Pierre Alliez, Florent Lafarge, and Andreas Fabri

#ifndef CGAL_URBAN_AREA_PROCESSING_INTERNAL_COMPARE_SCORES_H
#define CGAL_URBAN_AREA_PROCESSING_INTERNAL_COMPARE_SCORES_H

// #include <CGAL/license/Urban_area_processing.h>

// STL includes.
#include <vector>
#include <utility>
#include <algorithm>

// CGAL includes.
#include <CGAL/assertions.h>

// TODO:
// Remove this and substitute by the one from the Shape_detection package.

namespace CGAL {
namespace Urban_area_processing {
namespace internal {

  template<typename FT>
  struct Compare_scores {
    const std::vector<FT>& m_scores;
      
    Compare_scores(
      const std::vector<FT>& scores) : 
    m_scores(scores) 
    { }

    bool operator()(const std::size_t i, const std::size_t j) const {
      CGAL_assertion(i >= 0 && i < m_scores.size());
      CGAL_assertion(j >= 0 && j < m_scores.size());
      return m_scores[i] > m_scores[j];
    }
  };

} // internal
} // Urban_area_processing
} // CGAL

#endif // CGAL_URBAN_AREA_PROCESSING_INTERNAL_COMPARE_SCORES_H
