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

#ifndef CGAL_LEVELS_OF_DETAIL_INTERNAL_IMAGE_NEIGHBOR_QUERY_H
#define CGAL_LEVELS_OF_DETAIL_INTERNAL_IMAGE_NEIGHBOR_QUERY_H

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
  class Image_neighbor_query {

  public:
    using Traits = GeomTraits;
    using Size_pair = std::pair<std::size_t, std::size_t>;
    using Image = std::vector<Pixel>;
    using Idx_map = std::map<Size_pair, std::size_t>;

    Image_neighbor_query(
      const Image& image,
      const Idx_map& idx_map,
      const std::size_t num_labels,
      const bool version_4) :
    m_image(image),
    m_idx_map(idx_map),
    m_num_labels(num_labels),
    m_version_4(version_4)
    { }

    void operator()(
      const std::size_t query_index,
      std::vector<std::size_t>& neighbors) const {

      CGAL_assertion(query_index >= 0);
      CGAL_assertion(query_index < m_image.size());

      neighbors.clear();
      const Pixel& pixel = m_image[query_index];
      const std::size_t i = pixel.i;
      const std::size_t j = pixel.j;

      std::vector<std::size_t> ni, nj;
      if (m_version_4)
        get_grid_neighbors_4(i, j, ni, nj);
      else 
        get_grid_neighbors_8(i, j, ni, nj);

      for (std::size_t k = 0; k < ni.size(); ++k) {
        const std::size_t neighbor = m_idx_map.at(std::make_pair(ni[k], nj[k]));
        if (!m_image[neighbor].is_boundary)
          neighbors.push_back(neighbor);
      }
    }

  private:
    const Image& m_image;
    const Idx_map& m_idx_map;
    const std::size_t m_num_labels;
    const bool m_version_4;

    void get_grid_neighbors_4(
      const std::size_t i, const std::size_t j,
      std::vector<std::size_t>& ni, 
      std::vector<std::size_t>& nj) const {

      ni.clear(); nj.clear();
      ni.resize(4); nj.resize(4);

      CGAL_assertion(i > 0 && j > 0);

      ni[0] = i - 1; nj[0] = j;
      ni[1] = i;     nj[1] = j + 1;
      ni[2] = i + 1; nj[2] = j;
      ni[3] = i;     nj[3] = j - 1;
    }

    void get_grid_neighbors_8(
      const std::size_t i, const std::size_t j,
      std::vector<std::size_t>& ni, 
      std::vector<std::size_t>& nj) const {

      ni.clear(); nj.clear();
      ni.resize(8); nj.resize(8);

      CGAL_assertion(i > 0 && j > 0);

      ni[0] = i - 1; nj[0] = j - 1;
      ni[1] = i - 1; nj[1] = j;
      ni[2] = i - 1; nj[2] = j + 1;
      ni[3] = i;     nj[3] = j + 1;
      ni[4] = i + 1; nj[4] = j + 1;
      ni[5] = i + 1; nj[5] = j;
      ni[6] = i + 1; nj[6] = j - 1;
      ni[7] = i;     nj[7] = j - 1;
    }
  };

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_INTERNAL_IMAGE_NEIGHBOR_QUERY_H
