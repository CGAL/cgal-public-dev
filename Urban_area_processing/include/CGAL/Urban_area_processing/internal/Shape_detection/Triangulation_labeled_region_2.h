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

#ifndef CGAL_URBAN_AREA_PROCESSING_INTERNAL_TRIANGULATION_LABELED_REGION_2_H
#define CGAL_URBAN_AREA_PROCESSING_INTERNAL_TRIANGULATION_LABELED_REGION_2_H

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
  class Triangulation_labeled_region_2 {
  public:
    using Traits = GeomTraits;
    using Face_handle = FaceHandle;

    Triangulation_labeled_region_2(
      const std::vector<Face_handle>& faces,
      const std::size_t ref_label,
      const std::size_t min_faces) :
    m_faces(faces),
    m_ref_label(ref_label),
    m_min_faces(min_faces)
    { }

    bool is_part_of_region(
      const std::size_t,
      const std::size_t face_index, 
      const std::vector<std::size_t>& region) const {

      const bool is_valid = check_region(region);
      if (!is_valid) return false;

      const auto face = *(m_faces.begin() + face_index);
      return face->info().label == m_ref_label;
    }

    inline bool is_valid_region(
      const std::vector<std::size_t>& region) const {

      const bool is_valid = check_region(region);
      return is_valid && (region.size() >= m_min_faces);
    }

    void update(
      const std::vector<std::size_t>&) { }

  private:
    const std::vector<Face_handle>& m_faces;
    const std::size_t m_ref_label;
    const std::size_t m_min_faces;

    bool check_region(
      const std::vector<std::size_t>& region) const {
      
      if (region.size() == 0)
        return false;

      if (region.size() == 1) {
        const std::size_t face_index = region[0];
        const auto face = *(m_faces.begin() + face_index);
        return face->info().label == m_ref_label;
      }
      return true;
    }
  };

} // internal
} // Urban_area_processing
} // CGAL

#endif // CGAL_URBAN_AREA_PROCESSING_INTERNAL_TRIANGULATION_LABELED_REGION_2_H
