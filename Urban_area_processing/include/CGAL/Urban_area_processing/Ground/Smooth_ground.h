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

#ifndef CGAL_URBAN_AREA_PROCESSING_SMOOTH_GROUND_H
#define CGAL_URBAN_AREA_PROCESSING_SMOOTH_GROUND_H

// #include <CGAL/license/Urban_area_processing.h>

namespace CGAL {
namespace Urban_area_processing {

  /*!
    \ingroup PkgUrbanAreaProcessingRefGround

    \brief reconstructs a smooth ground from the input point cloud.

    This class creates a smooth surface on top of Delaunay triangulation, 
    refines this triangulation with respect to the local density of input points,
    and returns the refined surface bounded either by a bounding box or a convex hull.

    \tparam GeomTraits 
    must be a model of `Kernel`.

    \tparam InputRange
    must be a model of `ConstRange` whose iterator type is `RandomAccessIterator`.

    \tparam PointMap 
    must be an `LvaluePropertyMap` whose key type is the value type of the input 
    range and value type is `GeomTraits::Point_3`.

    \cgalModels `Ground`
  */
  template<
  typename GeomTraits,
  typename InputRange,
  typename PointMap>
  class Smooth_ground {

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
      defined in `CGAL::Urban_area_processing::Ground_parameters`
    */
    Smooth_ground(
      const InputRange& input_range,
      const PointMap point_map,
      const Ground_parameters<typename GeomTraits::FT> parameters) : 
    m_input_range(input_range),
    m_point_map(point_map),
    m_parameters(parameters) { 

      CGAL_precondition(input_range.size() > 0);
    }

    /// @}

    /// \name Reconstruction
    /// @{

    /*!  
      implements `Ground::reconstruct`
    */
    void reconstruct() {

    }

    /*!  
      creates a smooth ground from the planar ground.

      \param planar_ground
      an instance of the `CGAL::Urban_area_processing::Planar_ground` that represents 
      a planar ground object.
    */
    void convert_from(
      CGAL::Urban_area_processing::Planar_ground& planar_ground) {

    }

    /// @}

    /// \name Urban objects
    /// @{

    /*!  
      implements `Ground::place`
    */
    template<typename UrbanObject>
    void place(const UrbanObject& urban_object) {

    }

    /*!  
      implements `Ground::update`
    */
    void update() {

    }

    /// @}

    /// \name Export
    /// @{

    /*!  
      implements `Ground::output_to_city_gml()`
    */
    void output_to_city_gml(std::ofstream& os) {

    }
    
    /// @}

  private:
    const Input_range& m_input_range;
    const Point_map m_point_map;
    const Ground_parameters<typename GeomTraits::FT> m_parameters;
    
  };

} // Urban_area_processing
} // CGAL

#endif // CGAL_URBAN_AREA_PROCESSING_SMOOTH_GROUND_H
