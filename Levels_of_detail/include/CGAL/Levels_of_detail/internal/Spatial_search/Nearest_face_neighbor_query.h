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

#ifndef CGAL_LEVELS_OF_DETAIL_INTERNAL_NEAREST_FACE_NEIGHBOR_QUERY_H
#define CGAL_LEVELS_OF_DETAIL_INTERNAL_NEAREST_FACE_NEIGHBOR_QUERY_H

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
  class Nearest_face_neighbor_query {

  public:
    using Traits = GeomTraits;
    using Point_3 = typename Traits::Point_3;
    using Polygon = std::vector<Point_3>;
    using Indices = std::vector<std::size_t>;

    Nearest_face_neighbor_query(
      const std::vector<Polygon>& faces) :
    m_faces(faces) {

			m_neighbors.clear();
      m_neighbors.reserve(m_faces.size());

			Indices neighbors;
			for (std::size_t i = 0; i < m_faces.size(); ++i) {
				find_neighbors(m_faces[i], i, neighbors);
				m_neighbors.push_back(neighbors);
			}
    }

    void operator()(
      const std::size_t query_index,
      std::vector<std::size_t>& neighbors) const {

      CGAL_assertion(query_index >= 0);
      CGAL_assertion(query_index < m_neighbors.size());
      neighbors = m_neighbors[query_index];
    }

  private:
    const std::vector<Polygon>& m_faces;
    std::vector<Indices> m_neighbors;

    void find_neighbors(
      const Polygon& face,
      const std::size_t curr_index,
      std::vector<std::size_t>& neighbors) {

      neighbors.clear();
			for (std::size_t i = 0; i < m_faces.size(); ++i)
				if (i != curr_index && share_edge(face, m_faces[i]))
					neighbors.push_back(i);
    }

		bool share_edge(
      const Polygon& f1,
      const Polygon& f2) const {

      CGAL_assertion(f1.size() > 0);
      CGAL_assertion(f2.size() > 0);

			for (std::size_t i = 0; i < f1.size(); ++i) {
				const std::size_t ip = (i + 1) % f1.size();
				for (std::size_t j = 0; j < f2.size(); ++j) {
					const std::size_t jp = (j + 1) % f2.size();
					if (internal::are_equal_edges_3(
            f1[i], f1[ip], f2[j], f2[jp]))
						return true;
				}
			}
			return false;
		}
  };

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_INTERNAL_NEAREST_FACE_NEIGHBOR_QUERY_H
