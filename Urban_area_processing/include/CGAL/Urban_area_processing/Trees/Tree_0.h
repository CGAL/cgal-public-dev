// Copyright (c) 2020 SARL GeometryFactory (France).
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

#ifndef CGAL_URBAN_AREA_PROCESSING_TREE_0_H
#define CGAL_URBAN_AREA_PROCESSING_TREE_0_H

// #include <CGAL/license/Urban_area_processing.h>

namespace CGAL {
namespace Urban_area_processing {

  /*!
    \ingroup PkgUrbanAreaProcessingRefTrees

    \brief reconstructs a tree with respect to the City GML standard LOD0.

    This class projects tree points onto the xy plane, estimates a center point 
    of the tree, estimates a radius of this tree, and creates a triangulated disc 
    that represents the tree footprint.

    \tparam GeomTraits 
    must be a model of `Kernel`.

    \tparam InputRange
    must be a model of `ConstRange` whose iterator type is `RandomAccessIterator`.

    \tparam PointMap 
    must be an `LvaluePropertyMap` whose key type is the value type of the input 
    range and value type is `GeomTraits::Point_3`.
  */
  template<
  typename GeomTraits,
  typename InputRange,
  typename PointMap>
  class Tree_0 {

  public:
    /// \cond SKIP_IN_MANUAL
    using Traits = GeomTraits;
    using Input_range = InputRange;
    using Point_map = PointMap;
    /// \endcond

    /// \name Initialization
    /// @{

    /*!
      \brief initializes all internal data structures.
      
      \param input_range
      an instance of `InputRange` with 3D points

      \param point_map
      an instance of `PointMap` that maps an item from `input_range` 
      to `GeomTraits::Point_3`

      \param parameters
      defined in `CGAL::Urban_area_processing::Tree_parameters`
    */
    Tree_0(
      const InputRange& input_range,
      const PointMap point_map,
      const Tree_parameters<typename GeomTraits::FT> parameters) : 
    m_input_range(input_range),
    m_point_map(point_map),
    m_parameters(parameters) { 

      CGAL_precondition(input_range.size() > 0);
    }

    /// @}

    /// \name Reconstruction
    /// @{

    /*!  
      reconstructs a tree.

      \param tree
      an instance of the `CGAL::Urban_area_processing::Tree` that represents 
      a tree urban object.
    */
    void reconstruct(
      CGAL::Urban_area_processing::Tree& tree) {

    }

    /// @}

  private:
    const Input_range& m_input_range;
    const Point_map m_point_map;
    const Tree_parameters<typename GeomTraits::FT> m_parameters;

  };

} // Urban_area_processing
} // CGAL

#endif // CGAL_URBAN_AREA_PROCESSING_TREE_0_H
