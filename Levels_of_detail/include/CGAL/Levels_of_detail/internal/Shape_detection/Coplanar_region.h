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

#ifndef CGAL_LEVELS_OF_DETAIL_INTERNAL_COPLANAR_REGION
#define CGAL_LEVELS_OF_DETAIL_INTERNAL_COPLANAR_REGION

// STL includes.
#include <vector>

// CGAL includes.
#include <CGAL/assertions.h>

// Internal includes.
#include <CGAL/Levels_of_detail/internal/utils.h>

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

  template<typename GeomTraits>
  class Coplanar_region {

  public:
    using Traits = GeomTraits;
    using Point_3 = typename Traits::Point_3;
    using Polygon = std::vector<Point_3>;

    Coplanar_region(
      const std::vector<Polygon>& faces) :
    m_faces(faces)
    { }

    bool is_already_visited(
      const std::size_t,
      const std::size_t query_index,
      const bool is_visited) const { return false; }

    bool is_part_of_region(
      const std::size_t,
      const std::size_t query_index,
      const std::vector<std::size_t>& region) const {

      CGAL_precondition(region.size() > 0);
      CGAL_precondition(query_index >= 0);
      CGAL_precondition(query_index < m_faces.size());
      return are_coplanar(query_index, region[0]);
    }

    bool is_valid_region(const std::vector<std::size_t>& region) const {
      return region.size() >= 1;
    }

    void update(const std::vector<std::size_t>&) {
      // skipped!
    }

  private:
    const std::vector<Polygon>& m_faces;

    bool are_coplanar(
      const std::size_t idx1,
      const std::size_t idx2) const {

      const Polygon& f1 = m_faces[idx1];
			const Polygon& f2 = m_faces[idx2];

      CGAL_assertion(f1.size() > 0);
      CGAL_assertion(f2.size() > 0);

			std::size_t count = 0;
			for (std::size_t i = 0; i < f1.size(); ++i) {
				const std::size_t ip = (i + 1) % f1.size();
				const std::size_t ipp = (i + 2) % f1.size();

				for (std::size_t j = 0; j < f2.size(); ++j)
					if (internal::are_coplanar_3(
            f1[i], f1[ip], f1[ipp], f2[j]))
              ++count;
			}
			return count == f1.size() * f2.size();
    }
  };

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_INTERNAL_COPLANAR_REGION
