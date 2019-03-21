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
// Author(s)     : Dmitry Anisimov, Simon Giraudot, Pierre Alliez, Florent Lafarge, and Andreas Fabri

#ifndef CGAL_LEVELS_OF_DETAIL_INTERNAL_PROPERTY_MAP_H
#define CGAL_LEVELS_OF_DETAIL_INTERNAL_PROPERTY_MAP_H

#include <CGAL/license/Levels_of_detail.h>

// CGAL includes.
#include <CGAL/property_map.h>

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

  template<
  typename Item_range, 
  typename Property_map>
  struct Item_property_map {

    using key_type = std::size_t;
    using value_type = typename Property_map::value_type;
    using reference = const value_type&;
    using category = boost::lvalue_property_map_tag;

    const Item_range& m_item_range;
    const Property_map& m_property_map;

    Item_property_map(
      const Item_range& item_range, 
      const Property_map& property_map) : 
    m_item_range(item_range),
    m_property_map(property_map) 
    { }

    reference operator[](key_type item_index) const { 
                
      CGAL_precondition(item_index >= 0);
      CGAL_precondition(item_index < m_item_range.size());

      const auto& key = *(m_item_range.begin() + item_index);
      return get(m_property_map, key);
    }

    friend inline reference get(
      const Item_property_map& item_map, 
      key_type key) { 
      return item_map[key];
    }
  };

  template<
  typename Point_map_3,
  typename Point_2>
  struct Point_2_from_point_3_property_map {

    using key_type = typename Point_map_3::key_type;
    using value_type = Point_2;
    using reference = const value_type&;
    using category = boost::lvalue_property_map_tag;

    const Point_map_3& m_point_map_3;

    Point_2_from_point_3_property_map(
      const Point_map_3& point_map_3) : 
    m_point_map_3(point_map_3)
    { }

    friend reference get(
      const Point_2_from_point_3_property_map& pmap, 
      const key_type& key) {
      const auto& point_3 = get(pmap.m_point_map_3, key);
      return reinterpret_cast<reference>(point_3);
    }
  };

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_INTERNAL_PROPERTY_MAP_H
