// Copyright (c) 2019 INRIA Sophia-Antipolis (France).
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
//
// Author(s)     : Dmitry Anisimov, Simon Giraudot, Pierre Alliez, Florent Lafarge, and Andreas Fabri
//

#ifndef CGAL_LEVELS_OF_DETAIL_INTERNAL_SPACIAL_SEARCH_SPHERE_NEIGHBOR_QUERY_H
#define CGAL_LEVELS_OF_DETAIL_INTERNAL_SPACIAL_SEARCH_SPHERE_NEIGHBOR_QUERY_H

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
#include <CGAL/Fuzzy_sphere.h>
#include <CGAL/Search_traits_2.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/Search_traits_adapter.h>

// Internal includes.
#include <CGAL/Levels_of_detail/internal/property_map.h>

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

  template<
  typename GeomTraits,
  typename InputRange,
  typename PointMap>
  class Sphere_neighbor_query {

  public:
    using Traits = GeomTraits;
    using Input_range = InputRange;
    using Point_map = PointMap;

    using Point = typename Point_map::value_type;
    using FT = typename Traits::FT;

    using Index_to_point_map =
    internal::Item_property_map<Input_range, Point_map>;

    using Search_base = typename std::conditional<
      std::is_same<typename Traits::Point_2, Point>::value,
      CGAL::Search_traits_2<Traits>,
      CGAL::Search_traits_3<Traits> >::type;

    using Search_traits =
    CGAL::Search_traits_adapter<std::size_t, Index_to_point_map, Search_base>;

    using Splitter =
    CGAL::Sliding_midpoint<Search_traits>;

    using Fuzzy_sphere
    = CGAL::Fuzzy_sphere<Search_traits>;

    using Tree
    = CGAL::Kd_tree<Search_traits, Splitter, CGAL::Tag_true>;

    Sphere_neighbor_query(
      const Input_range& input_range,
      const FT sphere_radius,
      const Point_map& point_map) :
    m_input_range(input_range),
    m_sphere_radius(sphere_radius),
    m_point_map(point_map),
    m_index_to_point_map(m_input_range, m_point_map),
    m_tree(
      boost::counting_iterator<std::size_t>(0),
      boost::counting_iterator<std::size_t>(m_input_range.size()),
      Splitter(),
      Search_traits(m_index_to_point_map)) {

      CGAL_precondition(m_input_range.size() > 0);
      CGAL_precondition(m_sphere_radius >= FT(0));

      m_tree.build();
    }

    void operator()(
      const std::size_t query_index,
      std::vector<std::size_t>& neighbors) const {

      CGAL_precondition(query_index >= 0);
      CGAL_precondition(query_index < m_input_range.size());

      const std::size_t sphere_center = query_index;

      const Fuzzy_sphere sphere(
        sphere_center,
        m_sphere_radius,
        FT(0),
        m_tree.traits());

      neighbors.clear();
      m_tree.search(std::back_inserter(neighbors), sphere);
    }

    void operator()(
      const Point& sphere_center,
      std::vector<std::size_t>& neighbors) const {

      const Fuzzy_sphere sphere(
        sphere_center,
        m_sphere_radius,
        FT(0),
        m_tree.traits());

      neighbors.clear();
      m_tree.search(std::back_inserter(neighbors), sphere);
    }

    const Index_to_point_map& point_map() const {
      return m_index_to_point_map;
    }

  private:
    const Input_range& m_input_range;

    const FT m_sphere_radius;

    const Point_map& m_point_map;
    const Index_to_point_map m_index_to_point_map;

    Tree m_tree;
  };

} // namespace internal
} // namespace Levels_of_detail
} // namespace CGAL

#endif // CGAL_LEVELS_OF_DETAIL_INTERNAL_SPACIAL_SEARCH_SPHERE_NEIGHBOR_QUERY_H
