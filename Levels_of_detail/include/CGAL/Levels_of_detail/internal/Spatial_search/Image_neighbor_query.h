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
#include <memory>

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
    using Pixels = std::vector<Pixel>;
    using Idx_map = std::map<Size_pair, std::size_t>;
    using Indices = std::vector<std::size_t>;

    Image_neighbor_query(
      const Pixels& pixels,
      const Idx_map& idx_map,
      const bool use_boundaries) :
    m_pixels(pixels),
    m_idx_map(idx_map),
    m_use_boundaries(use_boundaries),
    m_use_version_4(true),
    m_use_version_8(false)
    { }

    void use_version_4() {
      m_use_version_4 = true;
      m_use_version_8 = false;

      m_ni.clear(); m_nj.clear();
      m_ni.resize(4); m_nj.resize(4);
    }

    void use_version_8() {
      m_use_version_4 = false;
      m_use_version_8 = true;

      m_ni.clear(); m_nj.clear();
      m_ni.resize(8); m_nj.resize(8);
    }

    void use_seeds(const Indices& seeds) {
      m_seeds = seeds;
    }

    void operator()(
      const std::size_t query_index,
      std::vector<std::size_t>& neighbors) {

      CGAL_assertion(query_index >= 0);
      CGAL_assertion(query_index < m_pixels.size());

      neighbors.clear();
      const Pixel& pixel = m_pixels[query_index];

      const std::size_t i = pixel.i;
      const std::size_t j = pixel.j;

      if (m_use_version_4)
        get_grid_neighbors_4(i, j, m_ni, m_nj);

      if (m_use_version_8)
        get_grid_neighbors_8(i, j, m_ni, m_nj);

      CGAL_assertion(m_ni.size() == m_nj.size());
      for (std::size_t k = 0; k < m_ni.size(); ++k) {
        const std::size_t neighbor =
          m_idx_map.at(std::make_pair(m_ni[k], m_nj[k]));

        if (m_seeds.size() > 0) {
          if (m_seeds[neighbor] != std::size_t(-1))
            neighbors.push_back(neighbor);
        } else {
          if (!m_use_boundaries) {
            if (m_pixels[neighbor].label != std::size_t(-1))
              neighbors.push_back(neighbor);
          } else
            neighbors.push_back(neighbor);
        }
      }
    }

  private:
    const Pixels& m_pixels;
    const Idx_map& m_idx_map;
    const bool m_use_boundaries;

    bool m_use_version_4;
    bool m_use_version_8;

    Indices m_ni, m_nj;
    Indices m_seeds;

    void get_grid_neighbors_4(
      const std::size_t i, const std::size_t j,
      Indices& ni,
      Indices& nj) const {

      CGAL_assertion(i > 0 && j > 0);

      ni[0] = i - 1; nj[0] = j;
      ni[1] = i;     nj[1] = j + 1;
      ni[2] = i + 1; nj[2] = j;
      ni[3] = i;     nj[3] = j - 1;
    }

    void get_grid_neighbors_8(
      const std::size_t i, const std::size_t j,
      Indices& ni,
      Indices& nj) const {

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
