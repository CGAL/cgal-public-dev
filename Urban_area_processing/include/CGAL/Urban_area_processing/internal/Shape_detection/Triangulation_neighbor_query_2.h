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

#ifndef CGAL_URBAN_AREA_PROCESSING_INTERNAL_TRIANGULATION_NEIGHBOR_QUERY_2_H
#define CGAL_URBAN_AREA_PROCESSING_INTERNAL_TRIANGULATION_NEIGHBOR_QUERY_2_H

// #include <CGAL/license/Urban_area_processing.h>

// STL includes.
#include <vector>
#include <utility>
#include <algorithm>

// CGAL includes.
#include <CGAL/assertions.h>

namespace CGAL {
namespace Urban_area_processing {
namespace internal {

  template<
  typename GeomTraits,
  typename FaceHandle>
  class Triangulation_neighbor_query_2 {
  public:
    using Traits = GeomTraits;
    using Face_handle = FaceHandle;

    Triangulation_neighbor_query_2(
      const std::vector<Face_handle>& faces) :
    m_faces(faces)
    { }

    void operator()(
      const std::size_t face_index,
      std::vector<std::size_t>& neighbors) const {

      neighbors.clear();
      const auto face = *(m_faces.begin() + face_index);
      for (std::size_t k = 0; k < 3; ++k) {
        const auto neighbor = face->neighbor(k);
        neighbors.push_back(neighbor->info().index);
      }
    }

  private:
    const std::vector<Face_handle>& m_faces;
  };

} // internal
} // Urban_area_processing
} // CGAL

#endif // CGAL_URBAN_AREA_PROCESSING_INTERNAL_TRIANGULATION_NEIGHBOR_QUERY_2_H
