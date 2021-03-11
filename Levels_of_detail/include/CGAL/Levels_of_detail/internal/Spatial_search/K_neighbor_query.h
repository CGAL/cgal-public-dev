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

#ifndef CGAL_LEVELS_OF_DETAIL_INTERNAL_SPACIAL_SEARCH_K_NEIGHBOR_QUERY_H
#define CGAL_LEVELS_OF_DETAIL_INTERNAL_SPACIAL_SEARCH_K_NEIGHBOR_QUERY_H

#include <CGAL/license/Levels_of_detail.h>

// STL includes.
#include <typeinfo>
#include <type_traits>

// Boost includes.
#include <CGAL/boost/iterator/counting_iterator.hpp>

// CGAL includes.
#include <CGAL/Kd_tree.h>
#include <CGAL/Splitters.h>
#include <CGAL/assertions.h>
#include <CGAL/Search_traits_2.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/Search_traits_adapter.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>

// Internal includes.
#include <CGAL/Levels_of_detail/internal/property_map.h>

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

template<
  typename GeomTraits,
  typename InputRange,
  typename PointMap>
  class K_neighbor_query {

  public:
    using Traits = GeomTraits;
    using Input_range = InputRange;
    using Point_map = PointMap;

    using FT = typename Traits::FT;

    using Point = typename Point_map::value_type;

    using Index_to_point_map =
    internal::Item_property_map<Input_range, Point_map>;

    using Search_base = typename std::conditional<
      std::is_same<typename Traits::Point_2, Point>::value,
      CGAL::Search_traits_2<Traits>,
      CGAL::Search_traits_3<Traits> >::type;

    using Search_traits =
    CGAL::Search_traits_adapter<std::size_t, Index_to_point_map, Search_base>;

    using Distance =
    CGAL::Distance_adapter<
      std::size_t,
      Index_to_point_map,
      CGAL::Euclidean_distance<Search_base> >;

    using Splitter =
    CGAL::Sliding_midpoint<Search_traits>;

    using Search_tree =
    CGAL::Kd_tree<Search_traits, Splitter, CGAL::Tag_true>;

    using Neighbor_search =
    CGAL::Orthogonal_k_neighbor_search<
      Search_traits,
      Distance,
      Splitter,
      Search_tree>;

    using Tree =
    typename Neighbor_search::Tree;

    K_neighbor_query(
      const Input_range& input_range,
      const FT k,
      const Point_map& point_map) :
    m_input_range(input_range),
    m_number_of_neighbors(static_cast<std::size_t>(CGAL::to_double(k))),
    m_point_map(point_map),
    m_index_to_point_map(m_input_range, m_point_map),
    m_distance(m_index_to_point_map),
    m_tree(
      boost::counting_iterator<std::size_t>(0),
      boost::counting_iterator<std::size_t>(m_input_range.size()),
      Splitter(),
      Search_traits(m_index_to_point_map)) {

      CGAL_precondition(m_input_range.size() > 0);
      CGAL_precondition(m_number_of_neighbors > 0);

      m_tree.build();
    }

    void operator()(
      const std::size_t query_index,
      std::vector<std::size_t>& neighbors) const {

      CGAL_precondition(query_index >= 0);
      CGAL_precondition(query_index < m_input_range.size());

      Neighbor_search neighbor_search(
        m_tree,
        get(m_index_to_point_map, query_index),
        m_number_of_neighbors,
        0,
        true,
        m_distance);

      neighbors.clear();
      for (auto it = neighbor_search.begin(); it != neighbor_search.end(); ++it)
        neighbors.push_back(it->first);
    }

    void operator()(
      const Point& query_point,
      std::vector<std::size_t>& neighbors) const {

      Neighbor_search neighbor_search(
        m_tree,
        query_point,
        m_number_of_neighbors,
        0,
        true,
        m_distance);

      neighbors.clear();
      for (auto it = neighbor_search.begin(); it != neighbor_search.end(); ++it)
        neighbors.push_back(it->first);
    }

    const Index_to_point_map& point_map() const {
      return m_index_to_point_map;
    }

  private:
    const Input_range& m_input_range;

    const std::size_t m_number_of_neighbors;

    const Point_map& m_point_map;
    const Index_to_point_map m_index_to_point_map;

    Distance m_distance;
    Tree m_tree;
  };

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_INTERNAL_SPACIAL_SEARCH_K_NEIGHBOR_QUERY_H
