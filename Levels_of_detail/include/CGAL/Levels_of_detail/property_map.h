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

#ifndef CGAL_LEVELS_OF_DETAIL_PROPERTY_MAP_H
#define CGAL_LEVELS_OF_DETAIL_PROPERTY_MAP_H

#include <CGAL/license/Levels_of_detail.h>

// CGAL includes.
#include <CGAL/property_map.h>

// Internal includes.
#include <CGAL/Levels_of_detail/enum.h>

namespace CGAL {
namespace Levels_of_detail {

  /*!
    \ingroup PkgLevelsOfDetailRef
      
    \brief Different property maps used by the `CGAL::Levels_of_detail::Levels_of_detail`.
  */

  /// \name Visibility
  /// @{

  /*!
    \ingroup PkgLevelsOfDetailRef

    \brief Maps a point to a visibility value in the range [0,1].

    \tparam SemanticMap 
    must be an `LvaluePropertyMap` whose key type is the value type of the 
    input range and value type is `CGAL::Levels_of_detail::Semantic_label`.
  */    
  template<typename SemanticMap>
  struct Visibility_from_semantic_map {

  public:
    
    /// \name Types
    /// @{
  
    /// Key type.
    typedef typename boost::property_traits<SemanticMap>::key_type key_type;
    
    /// Value type.
    typedef double value_type;

    /// @}

    /// \cond SKIP_IN_MANUAL
    using Semantic_map = SemanticMap;
    using reference = value_type;
    using category = boost::readable_property_map_tag;

    Semantic_map m_semantic_map;
    /// \endcond

    /// \name Constructors
    /// @{

    /*!
      \brief default constructor of a visibility map.
    */
    Visibility_from_semantic_map() { }

    /*!
      \brief initializes a visibility map with a semantic map.
    */
    Visibility_from_semantic_map(Semantic_map semantic_map) : 
    m_semantic_map(semantic_map) 
    { }

    /// @}

    /// \name Access
    /// @{

    /*!
      \brief returns a visibility value for a given key.
    */
    friend value_type get(
      const Visibility_from_semantic_map& visibility_map, 
      const key_type& key) {

      const Semantic_label label = get(visibility_map.m_semantic_map, key);
        
      if (label == Semantic_label::BUILDING_INTERIOR)
        return 1.0;
      if (label == Semantic_label::BUILDING_BOUNDARY)
        return 0.5;

      return 0.0; // ground, vegetation, unassigned classes
    }

    /// @}

  }; // Visibility_from_semantic_map

  /// @}

} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_PROPERTY_MAP_H
